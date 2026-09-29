"""General induction elimination, regular vorticity and advected exact-wave tests."""
import numpy as np
import pytest
from scipy.constants import atomic_mass,elementary_charge as E
from scipy.sparse import bmat
from scipy.sparse.linalg import splu
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.current_thermal import CurrentThermal
from plasma_bench.current_thermal_midpoint import CurrentThermalMidpoint
from plasma_bench.inertial_induction import InertialInduction
from plasma_bench.dual_kinematics import velocity
from plasma_bench.moving_current import MovingCurrent
from examples.moving_species_reference import MovingSpeciesReference
from examples.current_beltrami_reference import CurrentBeltrami


def test_general_magnetic_elimination_matches_original_mixed_system():
    ref=MovingSpeciesReference();mesh=MHDMesh(ref.radius,ref.length/2,8,16)
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');state=ref.fields(r,z,.2/ref.frequency)
    model=CurrentThermal(mesh,state['n'],ref.ion_mass);geo=model.geometry
    source=model.mass@geo.project_electric(ref.fields(geo.r,geo.z,.2/ref.frequency)['e'])
    actual=InertialInduction(model).solve(source)
    old=CurrentThermalMidpoint(model);tau=1e-8;impedance=model.inertia_coefficient.max()/tau;d=old.diagonal
    matrix=bmat([[old.scaled_mass,tau*impedance*d@old.stiffness@d],
                 [-d@model.inertia@d/(tau*impedance),old.scaled_mass]],format='csc')
    rhs=np.r_[np.zeros_like(source),old.scale*source/impedance];solved=splu(matrix).solve(rhs)
    electric=impedance*old.scale*solved[len(source):];direct=geo.faraday@electric
    error=actual.magnetic-direct
    assert np.sqrt(error@(model.hessian@error)/(direct@(model.hessian@direct)))<1e-10
    assert actual.residual<2e-12 and actual.iterations>0
    zero=InertialInduction(model).solve(np.zeros_like(source))
    np.testing.assert_array_equal(zero.magnetic,0.)


def test_regular_cylindrical_vorticity_includes_rotation_and_axis_terms():
    mesh=MHDMesh(.035,.12,5,10);n=np.full((6,10),1e20);model=CurrentThermal(mesh,n,39.948*atomic_mass)
    x=.5*mesh.rf[:,None]**2*np.ones((1,mesh.nz));regular=np.stack((2-4*x/mesh.radius**2,3+5*x,7+11*x),axis=-1)
    v,curl=velocity(model.geometry,regular);r=model.geometry.r;xx=.5*r*r
    expected=np.stack((r*(2-4*xx/mesh.radius**2),r*(3+5*xx),7+11*xx),axis=-1)
    np.testing.assert_allclose(v,expected,rtol=1e-14,atol=1e-15)
    np.testing.assert_allclose(curl,np.stack((np.zeros_like(r),-11*r,6+10*r*r),axis=-1),rtol=3e-11,atol=1e-12)


def test_galilean_boost_of_exact_inertial_wave_refines():
    ref=CurrentBeltrami(number_density=1e18,amplitude=1e-6);speed=5e4;errors=[]
    for nr in (4,8,16):
        mesh=MHDMesh(ref.radius,ref.length/2,nr,2*nr);n=np.full((nr+1,2*nr),ref.number_density)
        model=CurrentThermal(mesh,n,ref.ion_mass);psi,g=ref.initial_fields(mesh)
        ue=model.electron.from_pressure(n*ref.electron_ev*E);ui=model.ion.from_pressure(n*ref.ion_ev*E)
        v=np.zeros(n.shape+(3,));v[...,2]=speed
        rate=MovingCurrent(model).rate(psi,g,ue,ui,v,np.zeros_like(n))
        actual=model.geometry.magnetic_samples(*model.geometry.unpack(rate.induction.magnetic))
        exact=ref.fields(model.geometry.r,model.geometry.z)['magnetic_rate']*(1+ref.k*speed/ref.omega)
        weight=model.geometry.weight
        errors.append(np.sqrt(np.sum(weight[...,None]*(actual-exact)**2)/np.sum(weight[...,None]*exact**2)))
    assert errors[1]<.65*errors[0] and errors[2]<.65*errors[1],errors


def test_moving_rate_failure_is_atomic_and_rejects_incompatible_boundaries():
    from examples.moving_current_case import initialize
    from examples.moving_current_measurements import protocol
    _,model,psi,g,ue,ui,v,nt,source,_=initialize(protocol()[0])
    saved=[a.copy() for a in (psi,g,ue,ui,v,nt,source)]
    with pytest.raises(RuntimeError,match='did not converge'):
        MovingCurrent(model,max_iterations=1).rate(psi,g,ue,ui,v,nt,electric_source=source)
    for expected,actual in zip(saved,(psi,g,ue,ui,v,nt,source)):
        np.testing.assert_array_equal(expected,actual)
    vbad=v.copy();vbad[-1,:,0]=1.
    with pytest.raises(ValueError,match='Impermeable'):
        MovingCurrent(model).rate(psi,g,ue,ui,vbad,nt)
