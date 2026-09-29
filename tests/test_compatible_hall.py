"""Regular-axis, energy, implicit time and coupled Hall fluid checks."""
from dataclasses import replace
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0 as MU0, atomic_mass
from scipy.integrate import solve_ivp
from plasma_bench.mhd import MHDParameters
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.toroidal_state import from_primitive, primitive
from plasma_bench.toroidal_mhd import ToroidalMHD
from plasma_bench.compatible_magnetic import energy_density
from plasma_bench.hall_compatible import CompatibleHallCurl
from plasma_bench.hall_induction import hall_induction
from plasma_bench.hall_midpoint import HallMidpoint, HallConvergenceError
from plasma_bench.hall_fluid import HallMHD
from examples.hall_manufactured import HallManufactured
from examples.hall_work_snapshots import snapshot
from examples.compatible_hall_references import nonlinear_state, wave_state
from examples.compressible_hall_reference import CompressibleHallWave


def test_compatible_faraday_commutes_with_physical_curl():
    mesh = MHDMesh(.7, .6, 6, 8)
    geo = CompatibleHallCurl(mesh)
    electric = np.random.default_rng(554).normal(size=geo.offsets[-1])
    er, et, ez = geo.split_electric(electric)
    er = er.reshape(mesh.nr, mesh.nz)
    et = np.pad(et.reshape(mesh.nr-1, mesh.nz), ((1, 1), (0, 0)))
    ez = np.pad(ez.reshape(mesh.nr, mesh.nz), ((0, 1), (0, 0)))
    i, j, _, _ = np.indices(geo.r.shape)
    def field(r, z):
        x = (r*r-mesh.rf[i]**2)/(2*mesh.radial_measure[i])
        t = (z-mesh.zf[j])/mesh.dz
        a = (1-t)*et[i, j]+t*et[i, (j+1)%mesh.nz]
        b = (1-t)*et[i+1, j]+t*et[i+1, (j+1)%mesh.nz]
        return np.stack((r*((1-t)*er[i, j]+t*er[i, (j+1)%mesh.nz]),
                         ((1-x)*mesh.rf[i]*a+x*mesh.rf[i+1]*b)/r,
                         (1-x)*ez[i, j]+x*ez[i+1, j]), axis=-1)
    def differentiate(fn, x, h):
        return (fn(x-2*h)-8*fn(x-h)+8*fn(x+h)-fn(x+2*h))/(12*h)
    dr = differentiate(lambda r: field(r, geo.z), geo.r, 1e-5*mesh.dr)
    dz = differentiate(lambda z: field(geo.r, z), geo.z, 1e-5*mesh.dz)
    f = field(geo.r, geo.z)
    minus_curl = np.stack((dz[..., 1], dr[..., 2]-dz[..., 0], -dr[..., 1]-f[..., 1]/geo.r), axis=-1)
    dp, dg = geo.unpack(geo.faraday @ electric)
    actual = geo.magnetic_samples(dp, dg)
    assert np.linalg.norm(actual-minus_curl)/np.linalg.norm(actual) < 2e-10


def test_axis_and_global_nonlinear_hall_errors_converge_without_fitting():
    ref = HallManufactured()
    axis, means = [], []
    for n in (8, 16, 32):
        mesh = MHDMesh(ref.radius, ref.length/2, n, 2*n)
        geo = CompatibleHallCurl(mesh)
        exact = snapshot(ref, mesh)
        density = ref.fields(geo.r, geo.z)['number_density']
        out = hall_induction(geo, exact['psi'], exact['g'], density)
        axis.append(np.linalg.norm(out.g_rate[0]-exact['g_rate'][0])/np.linalg.norm(exact['g_rate'][0]))
        means.append(np.sqrt(np.sum(mesh.volumes*(out.g_rate-exact['g_rate'])**2)/np.sum(mesh.volumes*exact['g_rate']**2)))
        assert abs(out.magnetic_work_w)/out.absolute_work_scale_w < 2e-14
    assert axis[0]/axis[1] > 3.5 and axis[1]/axis[2] > 3.5
    assert means[0]/means[1] > 3.5 and means[1]/means[2] > 3.5
    assert axis[-1] < .002


def test_energy_hessian_and_quadrature_are_consistent_and_positive():
    ref = HallManufactured()
    mesh = MHDMesh(ref.radius, ref.length/2, 8, 12)
    geo = CompatibleHallCurl(mesh, 12)
    initial = snapshot(ref, mesh)
    stepper = HallMidpoint(geo)
    x = geo.pack(initial['psi'], initial['g'])
    derivative = geo.magnetic_gradient(initial['psi'], initial['g'])
    assert_allclose(stepper.hessian @ x, derivative, rtol=2e-12, atol=1e-11)
    actual = np.sum(mesh.volumes*geo.energy_density(initial['psi'], initial['g']))
    independent = np.sum(geo.weight*np.sum(geo.magnetic_samples(initial['psi'], initial['g'])**2, axis=-1))/(2*MU0)
    assert actual == pytest.approx(independent, rel=2e-14)
    assert .5*x @ (stepper.hessian @ x) == pytest.approx(actual, rel=3e-14)
    assert np.all(energy_density(mesh, initial['psi']) > 0)


@pytest.mark.parametrize('eta', [0., 1e-4])
def test_midpoint_pairs_magnetic_energy_poynting_and_positive_heat(eta):
    ref = HallManufactured()
    mesh = MHDMesh(ref.radius, ref.length/2, 8, 12)
    geo = CompatibleHallCurl(mesh)
    initial = snapshot(ref, mesh)
    density = ref.fields(geo.r, geo.z)['number_density']
    stepper = HallMidpoint(geo, tolerance=1e-12)
    out = stepper.step(initial['psi'], initial['g'], density, 2e-7, resistivity=eta)
    initial_energy = np.sum(mesh.volumes*geo.energy_density(initial['psi'], initial['g']))
    final_energy = np.sum(mesh.volumes*geo.energy_density(out.psi, out.g))
    assert abs(final_energy-initial_energy+out.physical_heat_j)/initial_energy < 2e-12
    assert abs(np.sum(mesh.volumes*out.total_energy_change))/initial_energy < 1e-15
    assert out.physical_heat_j >= 0
    assert (out.physical_heat_j > 0) == (eta > 0)
    assert out.nonlinear_residual < 1e-12
    assert out.linear_residual < 1e-12


def test_nonlinear_midpoint_time_order_against_independent_adaptive_rk():
    ref = replace(HallManufactured(), density_contrast=.2, poloidal_amplitude=.001)
    mesh = MHDMesh(ref.radius, ref.length/2, 4, 8)
    geo = CompatibleHallCurl(mesh)
    initial = snapshot(ref, mesh)
    density = ref.fields(geo.r, geo.z)['number_density']
    old = geo.pack(initial['psi'], initial['g'])
    scale = np.concatenate((np.full((mesh.nr+1)*mesh.nz, ref.radius**2), np.full(mesh.nr*mesh.nz, 1/ref.radius)))
    def rhs(_, y):
        psi, g = geo.unpack(y*scale)
        out = hall_induction(geo, psi, g, density)
        return geo.pack(out.psi_rate, out.g_rate)/scale
    end = 4e-7
    reference = solve_ivp(rhs, (0., end), old/scale, method='DOP853', rtol=2e-12, atol=2e-15)
    assert reference.success
    exact = reference.y[:, -1]*scale
    errors = []
    for count in (4, 8, 16):
        stepper = HallMidpoint(geo, tolerance=2e-12)
        psi, g = initial['psi'].copy(), initial['g'].copy()
        for _ in range(count):
            out = stepper.step(psi, g, density, end/count)
            psi, g = out.psi, out.g
        diff = geo.pack(psi, g)-exact
        errors.append(np.sqrt(diff @ (stepper.hessian @ diff)))
    assert 3.7 < errors[0]/errors[1] < 4.3
    assert 3.7 < errors[1]/errors[2] < 4.3


def test_midpoint_failure_keeps_input_state_unchanged():
    ref = HallManufactured()
    mesh = MHDMesh(ref.radius, ref.length/2, 6, 8)
    geo = CompatibleHallCurl(mesh)
    initial = snapshot(ref, mesh)
    psi, g = initial['psi'].copy(), initial['g'].copy()
    stepper = HallMidpoint(geo, tolerance=1e-13, max_iterations=1)
    with pytest.raises(HallConvergenceError):
        stepper.step(psi, g, ref.fields(geo.r, geo.z)['number_density'], 2e-7)
    assert_allclose(psi, initial['psi'], rtol=0., atol=0.)
    assert_allclose(g, initial['g'], rtol=0., atol=0.)


def test_midpoint_rejects_invalid_flux_before_zero_transport_shortcut():
    mesh = MHDMesh(1., 1., 4, 8)
    geo = CompatibleHallCurl(mesh)
    stepper = HallMidpoint(geo)
    r, _ = mesh.corners
    for where in ((0, 2), (-1, 2), (2, -1)):
        psi = .5*r*r
        psi[where] += .1
        with pytest.raises(ValueError):
            stepper.step(psi, np.zeros((4, 8)), None, .1)


def test_coupled_finite_pressure_fast_hall_wave_propagation_refines():
    wave = CompressibleHallWave(.035, .24, 1e20, 39.948*atomic_mass, .005, 16., amplitude=1.)
    errors = []
    for n in (6, 12, 24):
        mesh = MHDMesh(wave.radius, wave.length/2, n, 2*n)
        model = HallMHD(mesh, wave_state(mesh, wave), wave.ion_mass)
        end = 1/abs(wave.omega)
        model.advance(end, dt_max=end/24)
        expected = wave_state(mesh, wave, end)
        error = model.state.conserved[..., 5]-expected.conserved[..., 5]
        errors.append(np.sqrt(np.sum(mesh.volumes*error*error)/np.sum(mesh.volumes*expected.conserved[..., 5]**2)))
    assert errors[0]/errors[1] > 2.5 and errors[1]/errors[2] > 2.5
    assert errors[-1] < .01


def test_explicit_energy_selection_and_zero_hall_fluid_reduction():
    ref = HallManufactured()
    mesh = MHDMesh(ref.radius, ref.length/2, 6, 8)
    state = nonlinear_state(mesh, ref, 39.948*atomic_mass)
    assert_allclose(primitive(mesh, state, 5/3)[..., 3], 16., rtol=2e-15)
    legacy = from_primitive(mesh, state.psi, state.conserved[..., 0], 16., toroidal_over_radius=state.conserved[..., 5])
    with pytest.raises(ValueError):
        HallMHD(mesh, legacy, 39.948*atomic_mass)
    model = HallMHD(mesh, state, 39.948*atomic_mass, hall=False)
    a, _ = model.rhs(state)
    b, _ = ToroidalMHD(mesh, state).rhs(state)
    assert_allclose(a.conserved, b.conserved, rtol=0., atol=0.)
    assert_allclose(a.psi, b.psi, rtol=0., atol=0.)
    dt = .1*model.timestep_limit()
    reference = ToroidalMHD(mesh, state)
    reference.step(dt/2); reference.step(dt/2)
    model.step(dt)
    assert_allclose(model.state.conserved, reference.state.conserved, rtol=0., atol=0.)
    assert_allclose(model.state.psi, reference.state.psi, rtol=0., atol=0.)
    assert model.state.magnetic_quadrature == 'compatible'


def test_coupled_hall_updates_full_state_and_preserves_conservative_budgets():
    ref = HallManufactured()
    mesh = MHDMesh(ref.radius, ref.length/2, 6, 12)
    initial = nonlinear_state(mesh, ref, 39.948*atomic_mass)
    model = HallMHD(mesh, initial, 39.948*atomic_mass, MHDParameters(resistivity=1e-4))
    model.advance(3e-7, dt_max=1e-7)
    out = model.diagnostics()
    for key in ('mass_kg', 'total_energy_j'):
        assert abs(out[key+'_change']) < 2e-13*abs(out[key])
    assert out['min_pressure_pa'] > 0 and out['physical_ohmic_heat_j'] > 0
    assert np.linalg.norm(model.state.psi-initial.psi) > 0
    assert np.linalg.norm(model.state.conserved[..., 1:3]) > 0
    assert model.steps == 3 and out['hall_nonlinear_iterations'] > 3


def test_coupled_hall_failure_is_atomic(monkeypatch):
    ref = HallManufactured()
    mesh = MHDMesh(ref.radius, ref.length/2, 6, 8)
    initial = nonlinear_state(mesh, ref, 39.948*atomic_mass)
    model = HallMHD(mesh, initial, 39.948*atomic_mass)
    def fail(*args, **kwargs):
        raise HallConvergenceError('Injected nonlinear failure')
    monkeypatch.setattr(model.induction, 'step', fail)
    with pytest.raises(HallConvergenceError):
        model.step(1e-7)
    assert_allclose(model.state.conserved, initial.conserved, rtol=0., atol=0.)
    assert_allclose(model.state.psi, initial.psi, rtol=0., atol=0.)
    assert model.steps == model.time == model.physical_ohmic_heat == model.hall_iterations == 0


def test_continuum_wave_initializer_keeps_true_averages_and_energy_choice():
    wave = CompressibleHallWave(.035, .24, 1e20, 39.948*atomic_mass, .005, 16., amplitude=10.)
    mesh = MHDMesh(wave.radius, wave.length/2, 8, 16)
    low, high = wave_state(mesh, wave, order=10), wave_state(mesh, wave, order=14)
    assert_allclose(low.conserved, high.conserved, rtol=2e-13, atol=1e-22)
    assert low.magnetic_quadrature == 'compatible'
    assert np.min(primitive(mesh, low, wave.gamma)[..., 3]) > 0
