"""Moving resistive composition, actual source/time refinement and atomic failure."""
import numpy as np
import pytest
from scipy.constants import mu_0 as MU0
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.moving_plasma import MovingPlasma
from plasma_bench.material_resistivity import MaterialResistivity
from plasma_bench.resistive_moving import ResistiveMovingPlasma
import plasma_bench.resistive_moving as moving_module
from plasma_bench.waveforms import Gaussian


def test_zero_eta_matches_ideal_material_dynamics_exactly():
    d=MovingDomain(4,8);args=dict(initial_a=.001,pressure=.2,bias=.0003)
    ideal=MovingPlasma(d,**args);resistive=ResistiveMovingPlasma(d,resistivity=0.,**args)
    for _ in range(3):
        ideal.step(.001);resistive.step(.001)
    np.testing.assert_array_equal(ideal.state.coordinates,resistive.state.coordinates)
    np.testing.assert_array_equal(ideal.state.velocity,resistive.state.velocity)
    np.testing.assert_array_equal(ideal.flux.labels,resistive.flux.labels)
    assert resistive.diagnostics()["ohmic_heat_j"]==0


def test_nonideal_labels_inverse_is_exact_on_deformed_geometry():
    d=MovingDomain(4,8);s=MovingPlasma(d,0.,pressure=0.)
    x=d.fluid_reference.copy();x[:,0]*=1.1*(1+.1*x[:,1]);x[:,1]*=.9
    a=.001*(1+.1*x[:,0]+.2*x[:,1]);labels=s.flux.labels_for_coefficients(x,a)
    s.flux.labels=labels
    np.testing.assert_allclose(s.flux.coefficients(x),a,rtol=2e-15)


def test_resistive_moving_source_and_splitting_time_order():
    velocity=[];budgets=[]
    for dt in (.002,.001,.0005,.00025):
        d=MovingDomain(4,8);bias=Gaussian(np.sqrt(MU0),.02,.03);b=float(bias.value(0.))
        s=ResistiveMovingPlasma(d,lambda r,z:-.75*b*(1-r*r-z*z),
             pressure=lambda r,z:5.625*b*b/MU0*r*r*(1-r*r-z*z),
             bias=bias,resistivity=.1*MU0,max_timestep=dt)
        entropy=s.pressure_invariant.copy();mass=s.cell_mass.copy()
        for _ in range(round(.02/dt)):
            s.step(dt)
        velocity.append(s.state.velocity);budgets.append(abs(s.diagnostics()["response_energy_residual_j"]))
        assert np.all(s.pressure_invariant>entropy)
        assert np.array_equal(s.cell_mass,mass)
        assert s.cell_ohmic_heat.min()>0
        assert s.maximum_resistive_balance_residual<1e-11
    differences=[np.linalg.norm(velocity[i]-velocity[i+1]) for i in range(3)]
    assert np.all((np.array(differences[:-1])/differences[1:]>3.8)&(np.array(differences[:-1])/differences[1:]<4.2))
    assert np.all((np.array(budgets[:-1])/budgets[1:]>3.8)&(np.array(budgets[:-1])/budgets[1:]<4.2))


def test_second_resistive_half_step_failure_is_fully_atomic(monkeypatch):
    d=MovingDomain(4,8);s=ResistiveMovingPlasma(d,.001,pressure=.1,resistivity=.1*MU0)
    before=s.state.copy();labels=s.flux.labels.copy();entropy=s.pressure_invariant.copy();heat=s.cell_ohmic_heat.copy()
    original=moving_module.integrate_diffusion;calls=0
    def failing(operator,*args,**kwargs):
        nonlocal calls
        calls+=1
        if calls==2:
            raise ValueError("Injected second resistive failure")
        return original(operator,*args,**kwargs)
    monkeypatch.setattr(moving_module,"integrate_diffusion",failing)
    with pytest.raises(ValueError,match="second resistive"):
        s.step(.001)
    assert calls==2
    assert np.array_equal(s.state.coordinates,before.coordinates) and np.array_equal(s.state.velocity,before.velocity)
    assert np.array_equal(s.flux.labels,labels) and np.array_equal(s.pressure_invariant,entropy)
    assert np.array_equal(s.cell_ohmic_heat,heat)
    assert (s.time,s.steps,s.source_work,s.maximum_resistive_balance_residual)==(0.,0,0.,0.)
