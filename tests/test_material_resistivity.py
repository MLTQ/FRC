"""Continuum diffusion/current/heat convergence and discrete integration contracts."""
import numpy as np
import pytest
from scipy.constants import mu_0 as MU0
from scipy.linalg import eigh
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.moving_field import MovingField
from plasma_bench.material_resistivity import MaterialResistivity
from examples.spherical_decay_reference import SphericalDecay


def test_sphere_current_field_and_local_heat_converge_independently():
    ref=SphericalDecay();errors=[]
    for n in (4,8,16):
        d=MovingDomain(n,2*n);o=MaterialResistivity(d,d.fluid_reference,MU0*ref.diffusivity)
        a=ref.fields(*d.fluid_reference.T)["a"];heat=np.zeros(d.fluid_cell_count)
        initial=o.induced(a);energy0=np.pi/MU0*np.dot(initial,o.stiffness@initial)
        for _ in range(50):
            update=o.step(a,.002);a=update.material_a;heat+=update.cell_heat
            assert update.cell_heat.min()>=0
        e=o.elements[0];nc=d.fluid_cell_count;r,z=e.points[:nc,:,0],e.points[:nc,:,1]
        weight=2*np.pi*e.weights[:nc]*r;exact=ref.fields(r,z,.1)
        br,bz=MovingField.fields(e,update.induced)
        current=o.current_quadrature(o.current_over_radius(a))
        field=np.sqrt(np.sum(weight*((br[:nc]-exact["Br"])**2+(bz[:nc]-exact["Bz"])**2))/
                      np.sum(weight*(exact["Br"]**2+exact["Bz"]**2)))
        jerror=np.sqrt(np.sum(weight*(current-exact["Jtheta"])**2)/np.sum(weight*exact["Jtheta"]**2))
        target=np.sum(weight*ref.heat_density(r,z,.1),axis=1);volume=weight.sum(axis=1)
        qerror=np.sqrt(np.sum((heat-target)**2/volume)/np.sum(target**2/volume))
        errors.append((field,jerror,qerror))
        assert abs(update.final_energy+heat.sum()-energy0) < 1e-11
        assert update.vacuum_relative_residual < 1e-11
    errors=np.array(errors)
    assert np.all(errors[:-1,0]/errors[1:,0] > 1.9)
    assert np.all(errors[:-1,1:]/errors[1:,1:] > 3.8)


def test_time_order_against_exact_semidiscrete_mode():
    d=MovingDomain(4,8);o=MaterialResistivity(d,d.fluid_reference,.1*MU0)
    n=d.fluid_node_count;k=o.stiffness.toarray();m=o.mass.toarray()
    schur=k[:n,:n]-k[:n,n:]@np.linalg.solve(k[n:,n:],k[n:,:n])
    eigen,vectors=eigh(schur,m);initial=vectors[:,0]*np.sqrt(MU0)/vectors[0,0]
    exact=initial*np.exp(-.1*eigen[0]*.1);errors=[]
    for dt in (.01,.005,.0025):
        a=initial.copy()
        for _ in range(round(.1/dt)):
            a=o.step(a,dt).material_a
        errors.append(np.sqrt((a-exact)@m@(a-exact)))
    ratios=np.array(errors[:-1])/errors[1:]
    assert np.all((ratios>3.95)&(ratios<4.05))


def test_deformed_source_work_matches_independent_physical_energy():
    d=MovingDomain(5,10,quadrature=7);x=d.fluid_reference.copy()
    x[:,0]*=.95*(1+.04*x[:,1]);x[:,1]+=.02*x[:,0]**2
    o=MaterialResistivity(d,x,.1*MU0);a=.0003*(1-.2*x[:,0]**2+.1*x[:,1])
    u0=o.induced(a,.001);update=o.step(a,.02,bias_start=.001,bias_end=.0015)
    def physical(u):
        value=0.
        for e in o.elements:
            br,bz=MovingField.fields(e,u)
            value+=np.sum(e.weights*e.points[...,0]*(br*br+bz*bz))
        value+=np.dot(d.physical_exterior_diagonal,u[d.node_count:]**2)
        return np.pi/MU0*value
    residual=physical(update.induced)-physical(u0)+update.cell_heat.sum()-update.source_work
    assert abs(residual)<2e-10
    assert update.cell_heat.min()>0
    assert abs(update.balance_residual)<1e-11


def test_zero_resistivity_keeps_tiny_material_coefficients_and_has_no_heat():
    d=MovingDomain(4,8);o=MaterialResistivity(d,d.fluid_reference,0.)
    a=np.full(d.fluid_node_count,1e-22);u=o.step(a,.01,bias_start=.001,bias_end=.002)
    assert np.array_equal(a,u.material_a)
    assert np.all(u.cell_heat==0)
    assert abs(u.balance_residual)<1e-10


def test_invalid_resistive_inputs_are_rejected_without_mutating_input():
    d=MovingDomain(4,8)
    with pytest.raises(ValueError,match="resistivity"):
        MaterialResistivity(d,d.fluid_reference,-1.)
    o=MaterialResistivity(d,d.fluid_reference,.1*MU0);a=np.ones(d.fluid_node_count)*.001;before=a.copy()
    with pytest.raises(ValueError,match="positive"):
        o.step(a,-.1)
    assert np.array_equal(a,before)
