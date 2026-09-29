"""Independent nonlinear geometry/transport and toroidal diffusion verification."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd import AxisymmetricMHD, MHDParameters
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd_state import MHDState, NonphysicalState
from plasma_bench.mhd_benchmarks import radial_fast_mode
from plasma_bench.toroidal_state import ToroidalState, from_primitive, primitive, radius_squared
from plasma_bench.toroidal_mhd import ToroidalMHD
from plasma_bench.toroidal_induction import toroidal_resistivity
from examples.toroidal_references import equilibrium, rotating_contact, torsional_wave, toroidal_decay


def test_rigid_rotation_energy_uses_true_annular_moment():
    mesh = MHDMesh(1., 1., 8, 4)
    state = equilibrium(mesh)
    solver = ToroidalMHD(mesh, state)
    integral_r2 = np.pi*mesh.radius**4*mesh.half_length
    assert solver.diagnostics()['kinetic_energy_j'] == pytest.approx(.5*.4**2*integral_r2, rel=2e-14)
    assert solver.diagnostics()['angular_momentum_kg_m2_s'] == pytest.approx(.4*integral_r2, rel=2e-14)
    assert_allclose(primitive(mesh, state, 5/3)[..., 6], .4, rtol=2e-15)


@pytest.mark.parametrize('omega,current_slope', [(.7, 0.), (0., .5), (.4, .3)])
def test_exact_rotation_and_pinch_force_balance(omega, current_slope):
    mesh = MHDMesh(1., 1., 16, 4)
    state = equilibrium(mesh, omega=omega, current_slope=current_slope)
    solver = ToroidalMHD(mesh, state)
    rhs, _ = solver.rhs(state)
    assert np.max(abs(rhs.conserved)) < 2e-12
    solver.advance(.1)
    assert_allclose(solver.state.conserved, state.conserved, rtol=2e-12, atol=3e-14)


def test_unbalanced_force_matches_continuum_at_axis_and_wall():
    mesh = MHDMesh(1., 1., 24, 4)
    state = equilibrium(mesh, omega=.7, current_slope=.3, balanced=False)
    rhs, _ = ToroidalMHD(mesh, state).rhs(state)
    expected = (.7**2-2*.3**2)*mesh.radial_centroid[:, None]
    assert_allclose(rhs.conserved[..., 1], np.broadcast_to(expected, (24, 4)), atol=1e-13)


@pytest.mark.parametrize('order', [1, 2])
def test_zero_toroidal_components_recover_existing_poloidal_rhs(order):
    mesh = MHDMesh(1., 1., 12, 8)
    old, _ = radial_fast_mode(mesh, epsilon=.01, time=.04)
    new = ToroidalState(np.concatenate((old.conserved, np.zeros((12, 8, 2))), axis=-1), old.psi.copy())
    par = MHDParameters(spatial_order=order)
    a, _ = AxisymmetricMHD(mesh, old, par).rhs(old)
    b, _ = ToroidalMHD(mesh, new, par).rhs(new)
    assert_allclose(a.conserved, b.conserved[..., :4], rtol=1e-11, atol=1e-13)
    assert_allclose(a.psi, b.psi, rtol=1e-13, atol=1e-18)
    assert_allclose(b.conserved[..., 4:], 0., atol=0.)


def test_nonlinear_rotating_entropy_wave_transport_converges():
    errors = []
    for nz in (24, 48, 96):
        mesh = MHDMesh(1., 1., 8, nz)
        initial = rotating_contact(mesh)
        solver = ToroidalMHD(mesh, initial)
        solver.advance(.2)
        exact = rotating_contact(mesh, .2)
        fields = solver.state.conserved[..., [0, 4]]
        target = exact.conserved[..., [0, 4]]
        errors.append(np.sqrt(np.sum(mesh.volumes[..., None]*(fields-target)**2)/np.sum(mesh.volumes[..., None]*target**2)))
        out = solver.diagnostics()
        for key in ('mass_kg', 'angular_momentum_kg_m2_s', 'total_energy_j'):
            assert abs(out[key+'_change']) < 2e-13*abs(out[key])
        assert out['min_pressure_pa'] > 0
    assert errors[0]/errors[1] > 2.5
    assert errors[1]/errors[2] > 2.5
    assert errors[-1] < .001


def test_torsional_wave_converges_and_keeps_divergence():
    errors = []
    for nz in (24, 48, 96):
        mesh = MHDMesh(1., 1., 8, nz)
        solver = ToroidalMHD(mesh, torsional_wave(mesh))
        solver.advance(.2)
        exact = torsional_wave(mesh, .2)
        q, target = primitive(mesh, solver.state, 5/3), primitive(mesh, exact, 5/3)
        errors.append(np.linalg.norm(q[..., 6]-target[..., 6])/np.linalg.norm(target[..., 6]))
        divergence_scale = np.sqrt(MU0)*(1/mesh.dr+1/mesh.dz)
        assert solver.diagnostics()['max_face_divergence_t_m']/divergence_scale < 64*np.finfo(float).eps
    assert errors[0]/errors[1] > 2.5
    assert errors[1]/errors[2] > 2.5
    assert errors[-1] < .01


def test_toroidal_diffusion_has_positive_work_paired_heat():
    mesh = MHDMesh(1., 1., 17, 12)
    rng = np.random.default_rng(230)
    g = rng.normal(size=(mesh.nr, mesh.nz))*1e-4
    dg, sr, sz, heat = toroidal_resistivity(mesh, g, .13*MU0)
    loss = -np.sum(mesh.volumes*radius_squared(mesh)*g*dg)/MU0
    assert loss == pytest.approx(heat, rel=2e-14)
    assert heat > 0
    assert_allclose(sr[[0, -1]], 0., atol=0.)
    assert_allclose(sz[:, 0], sz[:, -1], atol=0.)
    # A finite regular current at the axis changes toroidal flux; it is not an exterior wall.
    a, b = -2e-4*mesh.radius**2, 1e-4
    g[:] = a+b*radius_squared(mesh)
    dg, _, _, _ = toroidal_resistivity(mesh, g, .13*MU0)
    assert_allclose(dg, 8*.13*b, rtol=3e-12)
    actual_flux_rate = np.sum(mesh.volumes*dg)/(2*np.pi)
    assert actual_flux_rate == pytest.approx(-2*.13*a*(2*mesh.half_length), rel=1e-13)


def test_toroidal_diffusion_matches_conducting_bessel_mode():
    errors = []
    for nr in (12, 24, 48):
        mesh = MHDMesh(1., 1., nr, 4)
        initial, _ = toroidal_decay(mesh)
        solver = ToroidalMHD(mesh, initial, MHDParameters(resistivity=.1*MU0))
        solver.advance(.1)
        exact, _ = toroidal_decay(mesh, .1)
        weight = mesh.volumes*radius_squared(mesh)
        errors.append(np.sqrt(np.sum(weight*(solver.state.conserved[..., 5]-exact.conserved[..., 5])**2)/np.sum(weight*exact.conserved[..., 5]**2)))
        out = solver.diagnostics()
        assert out['physical_ohmic_heat_j'] > 0
        assert abs(out['total_energy_j_change']) < 2e-13*out['total_energy_j']
    assert errors[0]/errors[1] > 2.5
    assert errors[1]/errors[2] > 2.5
    assert errors[-1] < .001


def test_rejected_trial_preserves_all_accepted_state(monkeypatch):
    import plasma_bench.toroidal_mhd as module
    mesh = MHDMesh(1., 1., 8, 8)
    solver = ToroidalMHD(mesh, rotating_contact(mesh))
    initial = solver.state.copy()
    original = module.primitive
    calls = 0
    def reject_final(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 5:
            raise NonphysicalState('Injected final-stage failure')
        return original(*args, **kwargs)
    monkeypatch.setattr(module, 'primitive', reject_final)
    with pytest.raises(NonphysicalState):
        solver.step()
    assert_allclose(solver.state.conserved, initial.conserved, rtol=0, atol=0)
    assert_allclose(solver.state.psi, initial.psi, rtol=0, atol=0)
    assert solver.time == solver.steps == solver.physical_ohmic_heat == 0


@pytest.mark.parametrize('kind', ['equilibrium', 'wave', 'decay'])
def test_reference_configuration_honors_nondefault_gamma(kind):
    from examples.toroidal_cases import ToroidalCase
    case = ToroidalCase(kind, gamma=1.4)
    solver = case.build()
    q = primitive(solver.mesh, solver.state, case.gamma)
    expected = 1.+(.5*case.omega**2-case.current_slope**2)*radius_squared(solver.mesh) if kind == 'equilibrium' else 1.
    assert_allclose(q[..., 3], np.broadcast_to(expected, q[..., 3].shape), rtol=3e-15)


def test_physical_field_error_includes_subcell_variation():
    from examples.toroidal_cases import ToroidalCase
    from examples.toroidal_field_errors import physical_errors
    errors = []
    for nz in (16, 32, 64):
        case = ToroidalCase('contact', nz=nz)
        solver = case.build()
        # Initial conservative averages are exact, but a limited reconstruction has spatial error.
        low, high = physical_errors(case, solver, 6), physical_errors(case, solver, 10)
        assert low['density'] == pytest.approx(high['density'], rel=1e-9)
        assert low['rotation_velocity'] == pytest.approx(high['rotation_velocity'], rel=1e-9)
        errors.append(high['density'])
        assert errors[-1] > 0
    assert errors[0]/errors[1] > 3
    assert errors[1]/errors[2] > 3
