"""Warm continuum species equations, coupled evolution and failure contracts."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.constants import atomic_mass, elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.warm_hall import WarmHallMHD
from plasma_bench.electron_midpoint import ElectronConvergenceError
from examples.warm_hall_reference import WarmHallWave
from examples.compressible_hall_reference import CompressibleHallWave
from examples.compatible_hall_references import wave_state


def reference(branch=2):
    return WarmHallWave(.035, .24, 1e20, 39.948*atomic_mass, .005, 6e20*E, branch=branch, amplitude=1.)


def initialize(ref, nr):
    mesh = MHDMesh(ref.radius, ref.length/2, nr, 2*nr)
    r, z = np.meshgrid(mesh.rf, mesh.zf[:-1], indexing='ij')
    pe = ref.electron_pressure+ref.fields(r, z)['electron_pressure']
    return WarmHallMHD(mesh, wave_state(mesh, ref), pe, ref.ion_mass, nonlinear_tolerance=1e-11)


@pytest.mark.parametrize('branch', [0, 1, 2])
def test_warm_wave_satisfies_both_species_momentum_thermal_and_maxwell(branch):
    ref = reference(branch)
    r, z = np.meshgrid(np.linspace(.002, .033, 11), np.linspace(-.11, .11, 13), indexing='ij')
    f = ref.fields(r, z, complex_output=True)
    h = ref.radius*1e-4
    def radial(key):
        samples = [ref.fields(r+i*h, z, complex_output=True)[key] for i in (-2, -1, 1, 2)]
        return (samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*h)
    def close(a, b):
        assert np.linalg.norm(a-b) < 4e-9*max(np.linalg.norm(a), np.linalg.norm(b))
    div = {}
    for key in ('velocity', 'electron_velocity'):
        v = f[key]
        div[key] = radial(key)[..., 0]+v[..., 0]/r+1j*ref.k*v[..., 2]
    close(-1j*ref.omega*f['density'], -ref.density*div['velocity'])
    close(-1j*ref.omega*f['electron_pressure'], -ref.gamma*ref.electron_pressure*div['electron_velocity'])
    close(-1j*ref.omega*f['ion_pressure'], -ref.gamma*(ref.pressure-ref.electron_pressure)*div['velocity'])
    bias = np.array([0., 0., ref.bias])
    pressure_gradient = f['electron_pressure_gradient']/ref.electron_fraction
    close(-1j*ref.omega*ref.density*f['velocity'], np.cross(f['current'], bias)-pressure_gradient)
    close(ref.number_density*E*(f['electric']+np.cross(f['electron_velocity'], bias)), -f['electron_pressure_gradient'])
    dr = radial('electric'); electric = f['electric']
    minus_curl = np.stack((1j*ref.k*electric[..., 1], dr[..., 2]-1j*ref.k*electric[..., 0],
                           -dr[..., 1]-electric[..., 1]/r), axis=-1)
    close(-1j*ref.omega*f['magnetic'], minus_curl)
    boundary = ref.fields(np.full_like(z[0], ref.radius), z[0], complex_output=True)
    assert np.max(abs(boundary['magnetic'][..., 1])) < 2e-14*np.max(abs(f['magnetic']))
    assert np.max(abs(boundary['electric'][..., 2])) > 0


def test_warm_reference_recovers_cold_limit_and_species_quadratic_energy():
    ref = reference()
    r, z = .017, .031
    fields = ref.fields(r, z, complex_output=True)
    species = abs(fields['electron_pressure'])**2/(ref.gamma*ref.electron_pressure)
    species += abs(fields['ion_pressure'])**2/(ref.gamma*(ref.pressure-ref.electron_pressure))
    assert species == pytest.approx(abs(fields['pressure'])**2/(ref.gamma*ref.pressure), rel=4e-16)
    cold = replace(ref, electron_fraction=0.)
    base = CompressibleHallWave(**{k: v for k, v in vars(cold).items() if k != 'electron_fraction' and k in CompressibleHallWave.__dataclass_fields__})
    np.testing.assert_array_equal(cold.fields(r,z)['electric'], base.fields(r,z)['electric'])


def test_actual_warm_hall_evolution_refines_and_preserves_energy():
    ref = reference(); errors = []
    for nr in (4, 8, 16):
        solver = initialize(ref, nr)
        solver.advance(.5/abs(ref.omega), dt_max=.025/abs(ref.omega))
        r, z = solver.thermal.mapping.geometry.r, solver.thermal.mapping.geometry.z
        exact = ref.fields(r, z, solver.time)['electron_pressure']
        actual = solver.thermal.pressure(solver.state, solver.electron_entropy)-ref.electron_pressure
        weight = solver.thermal.mapping.volumes
        errors.append(np.sqrt(np.sum(weight*(actual-exact)**2)/np.sum(weight*exact*exact)))
        diagnostics = solver.diagnostics()
        assert abs(diagnostics['total_energy_j_change']) < 3e-15*solver.initial_integrals['total_energy_j']
        assert diagnostics['min_electron_pressure_pa'] > 0 and diagnostics['min_ion_pressure_pa'] > 0
        assert diagnostics['physical_ohmic_heat_j'] == 0
    assert errors[1] < .35*errors[0] and errors[2] < .35*errors[1]


def test_failed_second_pressure_substep_preserves_every_accepted_account(monkeypatch):
    solver = initialize(reference(), 4)
    state, entropy, diagnostic = solver.state.copy(), solver.electron_entropy.copy(), solver.diagnostics()
    original, count = solver._pressure, 0
    def fail_second(*args):
        nonlocal count
        count += 1
        if count == 2:
            raise ElectronConvergenceError('forced second pressure failure')
        return original(*args)
    monkeypatch.setattr(solver, '_pressure', fail_second)
    with pytest.raises(ElectronConvergenceError, match='second pressure'):
        solver.step(.01/abs(reference().omega))
    np.testing.assert_array_equal(solver.state.conserved, state.conserved)
    np.testing.assert_array_equal(solver.state.psi, state.psi)
    np.testing.assert_array_equal(solver.electron_entropy, entropy)
    assert solver.diagnostics() == diagnostic
    with pytest.raises(ValueError, match='matching electron entropy'):
        solver.rhs(solver.state)
