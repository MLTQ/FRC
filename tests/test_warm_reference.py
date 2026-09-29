"""Continuum species/energy equations and actual finite-amplitude evolution."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.constants import mu_0 as MU0
from examples.warm_reference import WarmReference
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd import MHDParameters
from plasma_bench.toroidal_state import ToroidalState
from plasma_bench.warm_mhd import AdiabaticSpeciesMHD


def initialize(mesh, ref, *, cfl=.25, density_scale=1e20):
    means = ref.cell_means(mesh)
    values = np.zeros((mesh.nr, mesh.nz, 6))
    values[..., 0], values[..., 2], values[..., 3] = means['density'], means['momentum'], means['energy']
    r, _ = mesh.corners
    state = ToroidalState(values, ref.bias*r*r/2, 'compatible')
    pressure = np.broadcast_to(ref.fields(mesh.zf[:-1], 0)['electron_pressure'], (mesh.nr+1, mesh.nz))
    return AdiabaticSpeciesMHD(mesh, state, pressure, ref.ion_mass, MHDParameters(gamma=ref.gamma, cfl=cfl), density_scale=density_scale)


@pytest.mark.parametrize('kind', ['contact', 'simple_wave'])
@pytest.mark.parametrize('gamma', [5/3, 1.4])
def test_exact_references_satisfy_parent_species_and_total_equations(kind, gamma):
    ref = WarmReference(kind=kind, gamma=gamma)
    z = np.linspace(-.11, .11, 29)
    time = .37*ref.period
    dz, dt = 1e-4*ref.length, 1e-4*ref.period
    f = ref.fields(z, time)
    def derivative(name, *, temporal=False):
        samples = [ref.fields(z if temporal else z+i*dz, time+i*dt if temporal else time)[name] for i in (-2, -1, 1, 2)]
        return (samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*(dt if temporal else dz))
    dr, dv, dp = (derivative(n) for n in ('density', 'velocity', 'pressure'))
    residual = derivative('density', temporal=True)+f['velocity']*dr+f['density']*dv
    assert np.max(abs(residual))/(ref.rho0/ref.period) < 2e-10
    residual = derivative('velocity', temporal=True)+f['velocity']*dv+dp/f['density']
    assert np.max(abs(residual))/(ref.sound_speed/ref.period) < 3e-10
    for name in ('electron_pressure', 'ion_pressure'):
        residual = derivative(name, temporal=True)+f['velocity']*derivative(name)+gamma*f[name]*dv
        assert np.max(abs(residual))/(ref.p0/ref.period) < 3e-10
    flux_derivative = dv*(f['energy']+f['pressure']-ref.bias**2/(2*MU0))+f['velocity']*(derivative('energy')+dp)
    assert np.max(abs(derivative('energy', temporal=True)+flux_derivative))/(ref.p0/ref.period) < 5e-10


def test_exact_reference_quadrature_and_characteristic_domain():
    mesh = MHDMesh(.035, .12, 4, 17)
    ref = WarmReference(kind='simple_wave', wave_mach=.12)
    for k, a in ref.cell_means(mesh, .6*ref.period, 12).items():
        np.testing.assert_allclose(a, ref.cell_means(mesh, .6*ref.period, 18)[k], rtol=5e-15, atol=1e-19)
    with pytest.raises(ValueError, match='crossing'):
        ref.fields(np.array([0.]), 10*ref.period)


@pytest.mark.parametrize('kind', ['contact', 'simple_wave'])
def test_actual_species_evolution_refines_towards_exact_solution(kind):
    ref = WarmReference(kind=kind)
    errors = []
    for nz in (16, 32, 64):
        mesh = MHDMesh(.035, .12, 4, nz)
        solver = initialize(mesh, ref)
        solver.advance(.2*ref.period)
        pe, _ = solver.thermal.validate(solver.state, solver.electron_entropy)
        exact = ref.fields(mesh.zf[:-1], solver.time)['electron_pressure']
        errors.append(float(np.sqrt(np.mean((pe-exact)**2))/ref.pe0))
        assert abs(solver.diagnostics()['total_energy_j_change']) < 2e-16
    assert errors[1] < .45*errors[0]
    assert errors[2] < .45*errors[1]
