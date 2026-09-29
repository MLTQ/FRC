"""Continuum stationary and weak-wave controls under a magnetic background."""
import numpy as np
import pytest
from scipy.constants import atomic_mass, elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.current_thermal import CurrentThermal
from plasma_bench.current_thermal_midpoint import CurrentThermalMidpoint
from examples.current_beltrami_reference import CurrentBeltrami
from examples.current_thermal_norms import wave_norms
from tests.test_current_thermal import initialize


@pytest.mark.parametrize('nr,dt,bias', [(8,1e-8,.005),(8,1e-6,-.005),(16,1e-8,.005),(16,1e-6,.5)])
def test_uniform_continuum_equilibrium_survives_repeated_steps(nr,dt,bias):
    mesh=MHDMesh(.035,.12,nr,2*nr);n=np.full((nr+1,2*nr),1e20)
    model=CurrentThermal(mesh,n,39.948*atomic_mass);solver=CurrentThermalMidpoint(model)
    psi=.5*bias*mesh.corners[0]**2;g=np.zeros((nr,2*nr))
    ue=model.electron.from_pressure(n*5*E);ui=model.ion.from_pressure(n*E)
    initial=[a.copy() for a in (psi,g,ue,ui)]
    for _ in range(10):
        step=solver.step(psi,g,ue,ui,dt)
        psi,g,ue,ui=step.psi,step.g,step.electron_entropy,step.ion_entropy
        assert step.nonlinear_residual<1e-10 and step.stored_state_residual<1e-14
    # Independently differentiate stored flux: the exact physical answer is Bz=bias.
    bz=np.diff(psi,axis=0)/mesh.radial_measure[:,None]
    br=-np.diff(psi,axis=1)[1:]/(mesh.rf[1:,None]*mesh.dz)
    assert np.max(abs(bz-bias))<2e-13*abs(bias)
    assert max(np.max(abs(br)),np.max(abs(mesh.r[:,None]*g)))<2e-13*abs(bias)
    np.testing.assert_array_equal(ue,initial[2]);np.testing.assert_array_equal(ui,initial[3])


@pytest.mark.parametrize('amplitude', [1e-9,1e-12])
def test_tiny_physical_wave_converges_relative_to_perturbation_not_bias(amplitude):
    ref=CurrentBeltrami(number_density=1e18,amplitude=amplitude)
    errors=[]
    for nr in (4,8,16):
        model,psi,g,ue,ui=initialize(ref,nr);solver=CurrentThermalMidpoint(model)
        end=.5/abs(ref.omega)
        for _ in range(20):
            step=solver.step(psi,g,ue,ui,end/20)
            psi,g,ue,ui=step.psi,step.g,step.electron_entropy,step.ion_entropy
        metric=wave_norms(model.mesh,ref,end,psi,g,model.current(model.geometry.pack(psi,g)),order=8)
        errors.append(metric['magnetic_relative_L2'])
        assert metric['current_relative_L2']<.4
    assert errors[1]<.6*errors[0] and errors[2]<.6*errors[1]


def test_strong_background_audit_records_roundoff_and_rejects_changed_energy():
    from copy import deepcopy
    from examples.verify_current_precision import cases
    from examples.current_thermal_evolution import measure
    from examples.current_thermal_audit import endpoint
    case=next(c for c in cases() if c['name']=='weak_strong_bias')
    case.update(nr=4,nz=8,steps=5)
    record,archive=measure(case);mesh=MHDMesh(.035,.12,4,8)
    audit=endpoint(mesh,case,archive,'final',record['final'])
    assert 0<audit['pressure_replay_roundoff_bound_pa']<1e-6
    altered=deepcopy(archive);altered['final_total_energy']+=1.
    with pytest.raises(AssertionError):
        endpoint(mesh,case,altered,'final',record['final'])
