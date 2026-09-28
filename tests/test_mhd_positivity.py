import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd import AxisymmetricMHD, MHDParameters
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd_state import primitive, NonphysicalState
from plasma_bench.mhd_reconstruction import Reconstruction
from plasma_bench.mhd_equilibrium import PeriodicEquilibrium
from plasma_bench.mhd_benchmarks import radial_fast_mode, longitudinal_rarefaction, double_rarefaction_exact


def test_density_is_positive_at_multidimensional_corners():
    mesh = MHDMesh(1, 1, 12, 12)
    q = np.zeros((12, 12, 6)); q[..., 0] = 1; q[..., 3] = 1
    q[4, 5, 0] = q[5, 4, 0] = .01
    q[6, 5, 0] = q[5, 6, 0] = 3
    original = q.copy()
    with pytest.raises(NonphysicalState, match="Reconstructed"):
        Reconstruction(mesh, q, positivity=False)
    rec = Reconstruction(mesh, q)
    assert rec.statistics["density_limited_cells"] > 0
    for trace in (*rec.faces(), *rec.corners()):
        assert np.min(trace[..., 0]) > 0
    assert_array_equal(q, original)


def test_thermal_blend_and_geometric_source_are_consistent():
    mesh = MHDMesh(1, 1, 12, 24)
    initial = PeriodicEquilibrium(background_pressure=.1).state(mesh)
    original = initial.copy()
    rec = Reconstruction(mesh, primitive(mesh, initial, 5/3))
    assert rec.statistics["pressure_blended_cells"] > 0
    for trace in (*rec.faces(), *rec.corners()):
        assert np.min(trace[..., 3]) > 0
    # Independent higher-order quadrature of the actual primitive reconstruction
    # checks that the source integrates the modified p+B²/2mu, not stale pt.
    nodes, weights = np.polynomial.legendre.leggauss(20)
    integral = np.zeros((mesh.nr, mesh.nz))
    for node, weight in zip(nodes, weights):
        radius = mesh.r+.5*mesh.dr*node
        trace = rec.primitive(rec._raw_radial(radius), rec._auxiliary_pressure(radius))
        integral += weight*(trace[..., 3]+np.sum(trace[..., 4:]**2, axis=-1)/(2*MU0))
    integral *= .5*mesh.dr/mesh.radial_measure[:, None]
    assert_allclose(rec.radial_source(), integral, rtol=2e-14, atol=2e-13)
    assert_array_equal(initial.conserved, original.conserved)
    assert_array_equal(initial.psi, original.psi)


def test_limiter_leaves_resolved_smooth_solution_unchanged():
    mesh = MHDMesh(1, 1, 24, 4)
    initial, _ = radial_fast_mode(mesh)
    enabled = AxisymmetricMHD(mesh, initial)
    disabled = AxisymmetricMHD(mesh, initial, MHDParameters(positivity_reconstruction=False))
    for solver in (enabled, disabled):
        solver.advance(.04)
    assert_array_equal(enabled.state.conserved, disabled.state.conserved)
    assert_array_equal(enabled.state.psi, disabled.state.psi)
    assert enabled.diagnostics()["reconstruction_limiter"]["pressure_blended_cells"] == 0


def test_low_pressure_equilibrium_now_runs_conservatively_and_converges():
    flows, flux_errors, fractions = [], [], []
    for n in (12, 24, 48):
        mesh = MHDMesh(1, 1, n, 2*n)
        initial = PeriodicEquilibrium(background_pressure=.1).state(mesh)
        solver = AxisymmetricMHD(mesh, initial)
        before = dict(solver.limiter_totals)
        solver.rhs(initial)
        assert solver.limiter_totals == before  # Diagnostic evaluation is pure.
        solver.advance(.1)
        diag = solver.diagnostics()
        assert diag["min_pressure_pa"] > 0
        assert diag["min_density_kg_m3"] > 0
        assert abs(diag["mass_relative_change"]) < 1e-12
        assert abs(diag["total_energy_relative_change"]) < 1e-12
        field = max(np.max(abs(b)) for b in mesh.face_fields(solver.state.psi))
        assert diag["max_face_divergence_t_m"]*mesh.dr/field < 1e-13
        flows.append(diag["max_velocity_m_s"])
        flux_errors.append(np.linalg.norm(solver.state.psi-initial.psi)/np.linalg.norm(initial.psi))
        limiting = diag["reconstruction_limiter"]
        fractions.append(limiting["pressure_blended_cells"]/limiting["accepted_cell_stages"])
    assert flows[1] < flows[0] and flows[2] < .5*flows[1]
    assert flux_errors[1] < .5*flux_errors[0] and flux_errors[2] < .5*flux_errors[1]
    assert fractions[0] > fractions[1] > fractions[2]
    assert fractions[2] == 0


@pytest.mark.parametrize("scale", [.01, 100.])
def test_trace_bounds_are_relative_not_an_absolute_pressure_floor(scale):
    mesh = MHDMesh(1, 1, 12, 24)
    q = primitive(mesh, PeriodicEquilibrium(background_pressure=.1).state(mesh), 5/3)
    base = Reconstruction(mesh, q)
    q[..., 3] *= scale**2
    q[..., 4:] *= scale
    changed = Reconstruction(mesh, q)
    assert_allclose(changed.pressure_weight, base.pressure_weight, atol=3e-13, rtol=3e-12)
    for actual, expected in zip(changed.faces(), base.faces()):
        assert np.all(actual[..., 3] > 0)
        assert_allclose(actual[..., 3]/scale**2, expected[..., 3], atol=3e-13, rtol=3e-11)


def test_positivity_switch_is_explicit_and_typed():
    with pytest.raises(ValueError, match="boolean"):
        MHDParameters(positivity_reconstruction=1)


def test_magnetized_double_rarefaction_against_independent_fans():
    errors = []
    for n in (80, 160, 320):
        mesh = MHDMesh(1, 1, 4, n)
        initial = longitudinal_rarefaction(mesh)
        solver = AxisymmetricMHD(mesh, initial, MHDParameters(gamma=1.4))
        solver.advance(.1)
        q = primitive(mesh, solver.state, 1.4)
        rho, _, pressure = double_rarefaction_exact(mesh.z, .1)
        select = abs(mesh.z) < .5
        errors.append((np.mean(abs(q[1, select, 0]-rho[select])), np.mean(abs(q[1, select, 3]-pressure[select]))))
        assert q[..., 0].min() > 0 and q[..., 3].min() > 0
        assert_allclose(solver.state.psi, initial.psi, atol=2e-17, rtol=0)
        assert abs(solver.diagnostics()["total_energy_relative_change"]) < 1e-12
        assert abs(solver.diagnostics()["mass_relative_change"]) < 1e-12
    assert np.all(np.asarray(errors[1]) < .6*np.asarray(errors[0]))
    assert np.all(np.asarray(errors[2]) < .65*np.asarray(errors[1]))
    # L1 convergence does not ensure an accurate central minimum on these grids;
    # the report retains its relative error instead of claiming full resolution.
