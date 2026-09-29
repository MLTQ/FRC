"""Independent force, vacuum, interface and integral checks for analytic FRC."""
import numpy as np
from numpy.testing import assert_allclose
from numpy.polynomial.legendre import leggauss
from scipy.constants import mu_0 as MU0
from examples.spherical_frc_reference import SphericalFRCReference
from plasma_bench.geometry import Grid
from plasma_bench.topology import analyze_flux


def test_continuum_force_balance_and_vacuum_current_from_independent_derivatives():
    ref=SphericalFRCReference(.7,.003)
    r=np.array([.1,.25,.4,.6,1.0])*ref.radius
    z=np.array([.1,-.3,.2,.9,-1.1])*ref.radius
    h=1e-5*ref.radius
    p,br,bz,pressure,current=ref.fields(r,z)
    rp,rm=ref.fields(r+h,z),ref.fields(r-h,z)
    zp,zm=ref.fields(r,z+h),ref.fields(r,z-h)
    assert_allclose((zp[1]-zm[1]-rp[2]+rm[2])/(2*h*MU0),current,atol=2e-5,rtol=2e-9)
    assert_allclose((rp[3]-rm[3])/(2*h),current*bz,atol=1e-7,rtol=2e-8)
    assert_allclose((zp[3]-zm[3])/(2*h),-current*br,atol=1e-7,rtol=2e-8)
    assert_allclose((rp[0]-rm[0])/(2*h*r),bz,atol=1e-11,rtol=2e-9)
    assert_allclose(-(zp[0]-zm[0])/(2*h*r),br,atol=1e-11,rtol=2e-9)
    assert np.all(pressure>=0) and np.all(pressure[3:]==0) and np.all(current[3:]==0)


def test_interface_magnetic_continuity_pressure_and_normal_traction():
    ref=SphericalFRCReference()
    theta=np.linspace(0,np.pi,33)
    r,z=ref.radius*np.sin(theta),ref.radius*np.cos(theta)
    eps=1e-7
    inner=ref.fields((1-eps)*r,(1-eps)*z)
    outer=ref.fields((1+eps)*r,(1+eps)*z)
    assert_allclose(np.array(inner[1:3]),np.array(outer[1:3]),atol=2e-6*abs(ref.bias),rtol=0)
    exact_br=-1.5*ref.bias*np.sin(theta)*np.cos(theta)
    exact_bz=1.5*ref.bias*np.sin(theta)**2
    assert_allclose(exact_br*np.sin(theta)+exact_bz*np.cos(theta),0,atol=1e-18)
    traction_inner=inner[3]+(inner[1]**2+inner[2]**2)/(2*MU0)
    traction_outer=(outer[1]**2+outer[2]**2)/(2*MU0)
    assert_allclose(traction_inner,traction_outer,atol=3e-6*ref.bias**2/MU0,rtol=0)
    axis=ref.fields(0.,0.)
    assert axis[0]==axis[1]==axis[3]==axis[4]==0
    assert axis[2]/ref.bias==-1.5
    o=ref.fields(ref.radius/np.sqrt(2),0.)
    assert_allclose(o[1:3],0,atol=1e-18)
    assert o[3]>0


def test_independent_full_space_quadrature_matches_closed_form_energy():
    ref=SphericalFRCReference(.7,.003)
    nodes,weights=leggauss(24)
    angle=nodes[None,:]
    s=.5*ref.radius*(nodes[:,None]+1)
    r,z=s*np.sqrt(1-angle*angle),s*angle
    _,br,bz,p,_=ref.fields(r,z)
    volume=2*np.pi*s*s*(.5*ref.radius*weights[:,None])*weights[None,:]
    inside_total=np.sum((br*br+bz*bz)/(2*MU0)*volume)
    inside_induced=np.sum((br*br+(bz-ref.bias)**2)/(2*MU0)*volume)
    pressure_integral=np.sum(p*volume)
    # Independent infinite radial substitution x=a/s maps the entire vacuum
    # to x in (0,1); integrate its actual field, rather than adding a fitted tail.
    x=.5*(nodes[:,None]+1);s=ref.radius/x
    _,br,bz,p,j=ref.fields(s*np.sqrt(1-angle*angle),s*angle)
    volume=2*np.pi*ref.radius**3/x**4*(.5*weights[:,None])*weights[None,:]
    outside=np.sum((br*br+(bz-ref.bias)**2)/(2*MU0)*volume)
    exact=ref.integrals()
    assert_allclose([inside_total,inside_induced+outside,outside,pressure_integral],
                    [exact['inside_total_magnetic_energy_j'],exact['induced_field_energy_j'],
                     exact['exterior_induced_field_energy_j'],exact['pressure_integral_j']],rtol=3e-14)
    assert not p.any() and not j.any()


def test_existing_topology_detector_refines_true_spherical_separatrix():
    ref=SphericalFRCReference()
    errors=[]
    for n in (24,48,96):
        grid=Grid(1.4,1.5,n+1,2*n+1)
        topology=analyze_flux(grid,ref.fields(*grid.mesh)[0])
        assert topology["closed_flux"] and len(topology["o_points"])==1
        assert len(topology["axis_nulls"])==2 and not topology["x_points"]
        o=topology["o_points"][0]
        assert_allclose([o["r_m"],o["z_m"]],[1/np.sqrt(2),0],atol=2e-5)
        assert_allclose(sorted(p["z_m"] for p in topology["axis_nulls"]),[-1.,1.],atol=4e-5)
        exact=ref.integrals()["separatrix_volume_m3"]
        errors.append(abs(topology["separatrix_volume_m3"]/exact-1))
        assert topology["largest_closed_volume_m3"] < topology["separatrix_volume_m3"]
    assert errors[1]<errors[0]/3.4 and errors[2]<errors[1]/3.4
