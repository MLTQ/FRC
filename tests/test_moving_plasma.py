"""Independent material expansion, ideal labels, equilibrium and step contracts."""
import numpy as np
import pytest
from scipy.constants import mu_0 as MU0
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.moving_plasma import MovingPlasma,MovingState
from plasma_bench.waveforms import Gaussian


def test_exact_free_expanding_gas_sphere_converges():
    errors=[]
    for n in (4,8,16):
        d=MovingDomain(n,2*n)
        s=MovingPlasma(d,0.,pressure=lambda r,z:1-r*r-z*z,max_timestep=.001)
        s.advance(.02)
        scale=np.sqrt(1+2*.02**2)
        exact_velocity=2*.02/scale*d.fluid_reference
        error=np.sqrt(np.sum(s.node_mass[:,None]*(s.state.velocity-exact_velocity)**2)/
                      np.sum(s.node_mass[:,None]*exact_velocity**2))
        errors.append(error)
        assert s.diagnostics()["mass_kg"] == s.initial["mass_kg"]
        assert s.diagnostics()["maximum_boundary_radius_m"] > 1
        assert abs(s.diagnostics()["response_energy_residual_j"]) < 1e-7
        assert abs(s.diagnostics()["axial_momentum_kg_m_s"]) < 1e-13
    assert np.all(np.asarray(errors[:-1])/errors[1:] > 2.3)
    assert errors[-1] < .012


def test_spherical_frc_initial_force_error_converges():
    errors=[];b=np.sqrt(MU0)
    for n in (4,8,16):
        d=MovingDomain(n,2*n)
        s=MovingPlasma(d,lambda r,z:-.75*b*(1-r*r-z*z),
                       pressure=lambda r,z:5.625*r*r*(1-r*r-z*z),bias=b)
        acceleration=s.evaluate()["force"]/s.node_mass[:,None]
        errors.append(np.sqrt(np.sum(s.node_mass[:,None]*acceleration**2)/s.node_mass.sum()))
    assert np.all(np.asarray(errors[:-1])/errors[1:] > 2.7)


def test_ideal_flux_labels_and_axis_affine_limit():
    d=MovingDomain(6,12);b=np.sqrt(MU0)
    s=MovingPlasma(d,lambda r,z:.5*b*(1+.1*z),pressure=1.,bias=b)
    x=d.fluid_reference.copy();x[:,0]*=1.2;x[:,1]=.8*x[:,1]+.05
    a=s.flux.coefficients(x)
    np.testing.assert_allclose(a,s.flux.labels/1.2**2,rtol=1e-14)
    np.testing.assert_allclose(x[:,0]**2*a,d.fluid_reference[:,0]**2*s.flux.labels,rtol=1e-14)
    before=s.flux.labels.copy();s.step(.0001)
    assert np.array_equal(before,s.flux.labels)
    assert np.all(s.state.coordinates[d.axis,0] == 0)


def test_interior_material_flux_error_converges_for_nonuniform_deformation():
    errors=[]
    for n in (4,8,16):
        d=MovingDomain(n,2*n);x0=d.fluid_reference;x=x0.copy()
        x[:,0]*=1+.2*x0[:,1];x[:,1]+=.1*x0[:,0]**2
        a0=.5*(1-.2*x0[:,0]**2+.1*x0[:,1])
        s=MovingPlasma(d,a0,pressure=0.)
        a=s.flux.coefficients(x);c=d.cells[:d.fluid_cell_count]
        # Corresponding interior Q1 material points should carry the same psi.
        # This is deliberately not a test at the preserved vertices.
        r0=d.shape@x0[c,0].T;r=d.shape@x[c,0].T
        before=r0*r0*(d.shape@a0[c].T);after=r*r*(d.shape@a[c].T)
        e=d.elements(x0)[0];w=(e.weights[:len(c)]*e.points[:len(c),:,0]).T
        errors.append(np.sqrt(np.sum(w*(after-before)**2)/np.sum(w*before**2)))
    assert np.all(np.asarray(errors[:-1])/errors[1:] > 3.9)
    assert errors[-1] < .00024


def test_verlet_energy_error_has_second_order_time_refinement():
    errors=[]
    for dt in (.002,.001,.0005):
        d=MovingDomain(4,8)
        s=MovingPlasma(d,0.,pressure=lambda r,z:1-r*r-z*z)
        for _ in range(round(.04/dt)):
            s.step(dt)
        errors.append(abs(s.diagnostics()["response_energy_residual_j"]))
    assert np.all((np.asarray(errors[:-1])/errors[1:] > 3.8)&(np.asarray(errors[:-1])/errors[1:] < 4.2))


def test_rejected_moving_proposal_is_atomic():
    d=MovingDomain(4,8)
    s=MovingPlasma(d,0.,pressure=0.,velocity=10*d.fluid_reference)
    before=s.state.copy();mass=s.cell_mass.copy();labels=s.flux.labels.copy()
    with pytest.raises(ValueError,match="matching sphere"):
        s.step(1.)
    assert np.array_equal(s.state.coordinates,before.coordinates)
    assert np.array_equal(s.state.velocity,before.velocity)
    assert np.array_equal(s.cell_mass,mass) and np.array_equal(s.flux.labels,labels)
    assert (s.time,s.steps,s.source_work)==(0.,0,0.)


def test_source_requires_time_cap_and_material_pressure_is_not_repaired():
    d=MovingDomain(4,8)
    with pytest.raises(ValueError,match="source time cap"):
        MovingPlasma(d,0.,bias=Gaussian(.001,.1,.03))
    with pytest.raises(ValueError,match="no floors"):
        MovingPlasma(d,0.,pressure=-.1)
