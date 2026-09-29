"""Continuum current response, independent physical energy and comparison limits."""
import numpy as np
import pytest
from scipy.constants import mu_0 as MU0, electron_mass as ME, elementary_charge as E
from scipy.linalg import eigh
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.moving_field import MovingField
from plasma_bench.material_resistivity import MaterialResistivity
from plasma_bench.current_inertia import MaterialCurrent
from examples.inertial_sphere_reference import InertialSphere
from examples.inertial_screening_reference import InertialScreening
from examples.profiled_current_reference import ProfiledSphere


def test_inertial_reference_obeys_faraday_ampere_and_physical_ohm_law():
    ref=InertialSphere(electron_length=.4);r=np.array([.13,.41,.6]);z=np.array([.2,-.3,.4]);t=.13;h=1e-5
    a=ref.fields(r,z,t)
    curl=(ref.fields(r,z+h,t)['Br']-ref.fields(r,z-h,t)['Br']-ref.fields(r+h,z,t)['Bz']+ref.fields(r-h,z,t)['Bz'])/(2*h)
    np.testing.assert_allclose(curl/MU0,a['Jtheta'],rtol=2e-9)
    jdot=(ref.fields(r,z,t+h)['Jtheta']-ref.fields(r,z,t-h)['Jtheta'])/(2*h)
    adot=(ref.fields(r,z,t+h)['a']-ref.fields(r,z,t-h)['a'])/(2*h)
    np.testing.assert_allclose(-r*adot,a['Etheta'],rtol=2e-9)
    np.testing.assert_allclose(a['Etheta'],MU0*ref.diffusivity*a['Jtheta']+ME/(ref.electron_density*E**2)*jdot,rtol=2e-9)
    before=ref.integrals(0.);after=ref.integrals(t)
    assert after['total_response_energy_j']+after['ohmic_heat_j']==pytest.approx(before['total_response_energy_j'])


def test_actual_sphere_current_and_heat_converge_to_continuum():
    ref=InertialSphere();errors=[]
    for n in (4,8,16):
        d=MovingDomain(n,2*n);o=MaterialCurrent(d,d.fluid_reference,electron_density=ref.electron_density,resistivity=.1*MU0)
        a=ref.fields(*d.fluid_reference.T)['a'];heat=np.zeros(d.fluid_cell_count)
        for _ in range(25):
            step=o.step(a,.004);a=step.material_a;heat+=step.cell_heat
        e=o.elements[0];nc=d.fluid_cell_count;r,z=e.points[:nc,:,0],e.points[:nc,:,1];weight=2*np.pi*e.weights[:nc]*r
        exact=ref.fields(r,z,.1);current=o.current_quadrature(step.current_over_radius)
        jerr=np.sqrt(np.sum(weight*(current-exact['Jtheta'])**2)/np.sum(weight*exact['Jtheta']**2))
        target=np.sum(weight*ref.heat_density(r,z,.1),axis=1);volume=weight.sum(axis=1)
        qerr=np.sqrt(np.sum((heat-target)**2/volume)/np.sum(target**2/volume))
        errors.append((jerr,qerr));assert heat.min()>=0
        assert abs(step.balance_residual)<1e-18
    # Inertia retains the weak-current error in the interpolated initial flux:
    # surface derivative error limits current L2 to roughly h^1.5, not h^2.
    assert np.all(np.array(errors[:-1])/np.array(errors[1:])>2.7)


def test_collisionless_screening_continuum_and_low_carrier_limit():
    errors=[]
    for n in (4,8,16):
        ref=InertialScreening(electron_length=.4);d=MovingDomain(n,2*n)
        o=MaterialCurrent(d,d.fluid_reference,electron_density=ref.electron_density,resistivity=0.)
        step=o.step(np.zeros(d.fluid_node_count),.01,bias_end=ref.bias)
        e=o.elements[0];nc=d.fluid_cell_count;r,z=e.points[:nc,:,0],e.points[:nc,:,1];weight=2*np.pi*e.weights[:nc]*r
        exact=ref.fields(r,z);current=o.current_quadrature(step.current_over_radius)
        errors.append(np.sqrt(np.sum(weight*(current-exact['Jtheta'])**2)/np.sum(weight*exact['Jtheta']**2)))
        assert np.all(step.cell_heat==0)
        assert abs(step.final_magnetic_energy+step.final_current_energy-step.source_work)<1e-18
    assert np.all(np.array(errors[:-1])/errors[1:]>3.5)
    dilute=InertialScreening(electron_length=100.);r=.6;z=.2
    approximate=-dilute.electron_density*E**2*r*dilute.bias/(2*ME)
    assert dilute.fields(r,z)['Jtheta']==pytest.approx(approximate,rel=2e-5)


def test_variable_density_resistance_source_work_matches_physical_integrals():
    d=MovingDomain(5,10,quadrature=7);x=d.fluid_reference.copy();x[:,0]*=.95*(1+.02*x[:,1])
    ne=lambda r,z:1e15*(1+.5*r*r+.2*z)
    eta=lambda r,z:MU0*.1*(1+.3*z+.1*r*r)
    o=MaterialCurrent(d,x,electron_density=ne,resistivity=eta)
    a=1e-8*(1-.3*x[:,0]**2+.1*x[:,1]);step=o.step(a,.03,bias_start=1e-8,bias_end=1.5e-8)
    def physical(a,bias):
        u=o.induced(a,bias);q=o.current_over_radius(a,bias);value=0.
        for e in o.elements:
            br,bz=MovingField.fields(e,u);value+=np.sum(e.weights*e.points[...,0]*(br*br+bz*bz))*np.pi/MU0
        value+=np.pi/MU0*np.dot(d.physical_exterior_diagonal,u[d.node_count:]**2)
        e=o.elements[0];nc=d.fluid_cell_count;r,z=e.points[:nc,:,0],e.points[:nc,:,1]
        value+=np.sum(2*np.pi*r*e.weights[:nc]*ME*o.current_quadrature(q)**2/(2*E**2*ne(r,z)))
        return value
    residual=physical(step.material_a,1.5e-8)-physical(a,1e-8)+step.cell_heat.sum()-step.source_work
    assert abs(residual)<1e-18 and step.cell_heat.min()>0
    e=o.elements[0];nc=d.fluid_cell_count;r,z=e.points[:nc,:,0],e.points[:nc,:,1]
    expected=.03*np.sum(2*np.pi*r*e.weights[:nc]*eta(r,z)*o.current_quadrature(step.midpoint_current_over_radius)**2,axis=1)
    np.testing.assert_allclose(step.cell_heat,expected,rtol=1e-13)


def test_variable_profiles_converge_to_independent_continuum_mode():
    ref=ProfiledSphere();errors=[]
    for n in (4,8,16):
        d=MovingDomain(n,2*n);o=MaterialCurrent(d,d.fluid_reference,electron_density=ref.density,resistivity=ref.resistivity)
        a=ref.fields(*d.fluid_reference.T)['a'];heat=np.zeros(d.fluid_cell_count)
        for _ in range(20):
            step=o.step(a,.005);a=step.material_a;heat+=step.cell_heat
        e=o.elements[0];nc=d.fluid_cell_count;r,z=e.points[:nc,:,0],e.points[:nc,:,1];w=2*np.pi*e.weights[:nc]*r
        exact=ref.fields(r,z,.1)['Jtheta'];j=o.current_quadrature(step.current_over_radius)
        qtarget=np.sum(w*ref.heat_density(r,z,.1),axis=1);vol=w.sum(axis=1)
        errors.append((np.sqrt(np.sum(w*(j-exact)**2)/np.sum(w*exact**2)),np.sqrt(np.sum((heat-qtarget)**2/vol)/np.sum(qtarget*qtarget/vol))))
    assert np.all(np.array(errors[:-1])/np.array(errors[1:])>2.7)


def test_massless_comparison_recovers_existing_resistive_component():
    d=MovingDomain(4,8);a=1e-8*(1-.2*d.fluid_reference[:,0]**2)
    o=MaterialCurrent(d,d.fluid_reference,electron_density=None,resistivity=.1*MU0)
    old=MaterialResistivity(d,d.fluid_reference,.1*MU0)
    x=o.step(a,.02,bias_start=1e-8,bias_end=2e-8);y=old.step(a,.02,bias_start=1e-8,bias_end=2e-8)
    np.testing.assert_allclose(x.material_a,y.material_a,rtol=1e-12,atol=1e-22)
    np.testing.assert_allclose(x.cell_heat,y.cell_heat,rtol=1e-11,atol=1e-14*y.cell_heat.max())
    assert x.final_current_energy==0


def test_current_time_order_uses_independently_diagonalized_mode():
    d=MovingDomain(4,8);o=MaterialCurrent(d,d.fluid_reference,electron_density=ME/(MU0*E**2*.2**2),resistivity=.1*MU0)
    n=d.fluid_node_count;k=o.stiffness.toarray();m=o.mass.toarray()
    s=k[:n,:n]-k[:n,n:]@np.linalg.solve(k[n:,n:],k[n:,:n]);e,v=eigh(s,m)
    initial=v[:,0]*1e-8;rate=.1*e[0]/(1+.2**2*e[0]);exact=initial*np.exp(-rate*.2);errors=[]
    for dt in (.04,.02,.01):
        a=initial.copy()
        for _ in range(round(.2/dt)):a=o.step(a,dt).material_a
        errors.append(np.sqrt((a-exact)@m@(a-exact)))
    np.testing.assert_allclose(np.array(errors[:-1])/errors[1:],4,rtol=.002)


def test_massless_ideal_limit_preserves_tiny_flux_under_large_source_change():
    d=MovingDomain(4,8);o=MaterialCurrent(d,d.fluid_reference,electron_density=None,resistivity=0.)
    a=np.full(d.fluid_node_count,1e-22)
    step=o.step(a,.01,bias_start=.001,bias_end=.002)
    assert np.array_equal(step.material_a,a)
    assert step.final_current_energy==0 and np.all(step.cell_heat==0)
    assert step.ampere_relative_residual<1e-12
    np.testing.assert_array_equal(step.induced,o.induced(a,.002))


def test_invalid_internal_vacuum_and_failed_step_do_not_mutate_fields():
    d=MovingDomain(4,8)
    for density in (0.,-1.,np.nan):
        with pytest.raises(ValueError):MaterialCurrent(d,d.fluid_reference,electron_density=density,resistivity=0.)
    o=MaterialCurrent(d,d.fluid_reference,electron_density=1e15,resistivity=0.);a=np.zeros(d.fluid_node_count)
    with pytest.raises(ValueError):o.step(a,-1.)
    assert np.array_equal(a,np.zeros_like(a))
    step=o.step(a,.1)
    assert np.array_equal(step.material_a,a) and step.balance_residual==0.
    stationary=np.full(d.fluid_node_count,1e-8);step=o.step(stationary,.1)
    assert np.array_equal(step.material_a,stationary)
