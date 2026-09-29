"""Physical scalar metric, thermodynamic chain rule and pressure-work verification."""
import numpy as np
import pytest
from scipy.constants import atomic_mass,elementary_charge as E,mu_0 as MU0
from scipy.sparse import kron,diags,csr_matrix
from scipy.sparse.linalg import spsolve
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.variable_compatible_current import VariableCompatibleCurrent
from plasma_bench.entropy_scalar import EntropyScalar
from plasma_bench.compatible_log_thermal import CompatibleLogThermal
from plasma_bench.kernel_induction import KernelInduction
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from examples.finite_mass_manufactured import cell_average
from examples.compatible_field_norms import direct_basis
from examples.variable_current_reference import VariableCurrentReference


def initialize(nr=8):
    ref=VariableCurrentReference();mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr)
    n=cell_average(mesh,lambda r,z:ref.fields(r,z)['n'],order=12)
    model=VariableCompatibleCurrent(mesh,n,ref.ion_mass)
    psi=ref.fields(*mesh.corners)['psi'];psi[[0,-1]]=0.;psi[:,-1]=psi[:,0]
    g=cell_average(mesh,lambda r,z:ref.fields(r,z)['g'],order=12)
    field=HarmonicMagnetic(mesh,psi,g,ref.axial_bias)
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij')
    return ref,model,field,r,z


def test_scalar_physical_polynomial_projection_and_independent_dense_metric():
    ref,m,field,r,z=initialize(6);geo=m.geometry;space=EntropyScalar(m.density_map)
    scalar=1+(.5*r*r)**2;physical=1+(.5*geo.r**2)**2
    np.testing.assert_allclose(space.sample(scalar),physical,rtol=2e-15)
    projected=space.project(physical)
    np.testing.assert_allclose(projected.values,scalar,rtol=2e-13,atol=2e-14)
    # Independent root-product basis and direct physical weighted Gram solve.
    nx,_,_=direct_basis(.5*m.mesh.rf**2,.5*geo.radial_points**2,5)
    nz,_,_=direct_basis(m.mesh.zf,geo.axial_points,5,True)
    N=kron(csr_matrix(nx),nz,format='csr')
    weights=geo._unlayout(geo.weight*m.density).ravel()
    mass=(N.T@diags(weights)@N).tocsc()
    rng=np.random.default_rng(72);rhs=rng.normal(size=space.shape)*1e10
    result=space.solve(rhs)
    target=spsolve(mass,rhs.ravel()).reshape(space.shape)
    np.testing.assert_allclose(result.values,target,rtol=1e-10,atol=1e-16)
    assert result.residual<space.tolerance


@pytest.mark.parametrize('gamma',(5/3,1.4))
def test_physical_entropy_and_pressure_work_with_original_induction(gamma):
    ref,m,field,r,z=initialize();thermal=CompatibleLogThermal(m,gamma=gamma)
    exact=ref.thermal(r,z,gamma=gamma)
    s=np.log(exact['entropy']/ref.fields(r,z)['n']/E)
    current=field.current(m);pair=thermal.pair(current,s)
    zero=np.zeros_like(m.density_map.cell_means)
    rate=thermal.energy_rate(s,pair.entropy_rate,zero)['total']
    work=current@(m.mass@pair.electric_full)
    measured=np.sum(m.mesh.volumes*rate)
    scale=abs(current)@abs(m.mass@pair.electric_full)+np.sum(m.mesh.volumes*abs(rate))
    assert abs(measured-work)<5e-12*scale
    physical_entropy=thermal.entropy_rate(s,pair.entropy_rate,zero)
    assert abs(np.sum(m.mesh.volumes*physical_entropy))<5e-12*np.sum(m.mesh.volumes*abs(physical_entropy))
    solved=KernelInduction(m,tolerance=1e-10,max_iterations=2000).solve(m.mass@pair.electric_reduced)
    kinetic=m.current_energy_rate(current,solved.current,zero)['total']
    magnetic=m.magnetic_energy_rate(field,solved.magnetic)
    assert abs(np.sum(m.mesh.volumes*(magnetic+kinetic+rate)))<1e-10*scale


def test_physical_density_and_entropy_energy_chain_rule_without_state_repair():
    ref,m,field,r,z=initialize();thermal=CompatibleLogThermal(m)
    s=1+.3*(r/m.mesh.radius)**2+.2*np.sin(np.pi*z/m.mesh.half_length)
    st=.4+.1*np.cos(np.pi*z/m.mesh.half_length)
    nt=.3*m.density_map.cell_means
    rate=thermal.energy_rate(s,st,nt)
    np.testing.assert_allclose(np.sum(m.mesh.volumes*rate['density']),
        np.sum(thermal.density_energy_gradient(s)*nt),rtol=2e-14)
    errors=[]
    for step in (.02,.01,.005):
        values=[]
        for sign in (-1,1):
            changed=VariableCompatibleCurrent(m.mesh,m.density_map.cell_means+sign*step*nt,m.ion_mass,geometry=m.geometry)
            other=CompatibleLogThermal(changed)
            values.append(other.energy_density(s+sign*step*st))
        finite=(values[1]-values[0])/(2*step)
        errors.append(np.linalg.norm(finite-rate['total'])/np.linalg.norm(rate['total']))
    assert errors[1]<.26*errors[0] and errors[2]<.26*errors[1]
    assert errors[-1]<3e-6


def test_constant_entropy_no_false_current_transport_and_atomic_extreme_rejection():
    ref,m,field,r,z=initialize();thermal=CompatibleLogThermal(m)
    s=np.full_like(r,2.);saved=s.copy();current=field.current(m)
    pair=thermal.pair(current,s)
    assert np.count_nonzero(pair.entropy_rate)==0
    assert np.count_nonzero(pair.electric_reduced)==0
    gradient_curl=m.geometry.faraday@pair.electric_full
    bound=abs(m.geometry.faraday)@abs(pair.electric_full)
    assert np.linalg.norm(gradient_curl)<1e-14*np.linalg.norm(bound)
    for value in (1000.,-1000.):
        with pytest.raises(ValueError,match='no floor'):thermal.pair(current,np.full_like(s,value))
    np.testing.assert_array_equal(s,saved)


def test_smooth_continuum_log_entropy_advection_refines_at_axis():
    errors=[];axis=[]
    for nr in (8,16,32):
        ref,m,field,r,z=initialize(nr);thermal=CompatibleLogThermal(m)
        exact=ref.thermal(r,z);s=np.log(exact['entropy']/ref.fields(r,z)['n']/E)
        pair=thermal.pair(field.current(m),s)
        wanted=exact['entropy_rate']/exact['entropy']
        errors.append(np.linalg.norm(pair.entropy_rate-wanted)/np.linalg.norm(wanted))
        axis.append(np.linalg.norm(pair.entropy_rate[:4]-wanted[:4])/np.linalg.norm(wanted[:4]))
    assert errors[1]<.4*errors[0] and errors[2]<.4*errors[1],errors
    assert axis[1]<.4*axis[0] and axis[2]<.4*axis[1],axis


def test_independent_polynomial_barotropic_pressure_force():
    mesh=MHDMesh(.7,.6,8,12);n0=1e18;nref=1e20;s0=.3
    n=cell_average(mesh,lambda r,z:n0*(1+(r/mesh.radius)**2))
    m=VariableCompatibleCurrent(mesh,n,39.948*atomic_mass)
    # gamma=2 is an explicit polynomial verification fixture, not an argon
    # closure. h=(n/nref)e exp(s0) and E=−2 grad(h)/e are exact polynomials.
    thermal=CompatibleLogThermal(m,gamma=2.,density_scale=nref)
    pair=thermal.pair(np.zeros(m.mass.shape[0]),np.full((mesh.nr+1,mesh.nz),s0))
    actual=m.geometry.sample_electric(pair.electric_full);r=m.geometry.r
    wanted=np.stack((-4*np.exp(s0)*n0/nref*r/mesh.radius**2,np.zeros_like(r),np.zeros_like(r)),axis=-1)
    np.testing.assert_allclose(actual,wanted,rtol=2e-10,atol=1e-12)
