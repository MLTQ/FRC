import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd_state import from_primitive, primitive, NonphysicalState
from plasma_bench.mhd import AxisymmetricMHD, MHDParameters
from plasma_bench.mhd_flux import ghost_cells, fast_speed
from plasma_bench.mhd_induction import ideal_emf
from plasma_bench.mhd_benchmarks import radial_fast_mode, advected_flux, resistive_bessel_mode, sod_initial, sod_exact, radial_pressure_balance, radial_average


def test_annular_volumes_and_discrete_field_divergence():
    mesh = MHDMesh(1, 2, 20, 40)
    assert_allclose(mesh.volumes.sum(), 4*np.pi)
    r, z = mesh.corners
    psi = r*r*(1-r*r)**3*np.cos(np.pi*z/2)
    bscale = max(np.abs(x).max() for x in mesh.face_fields(psi))
    assert np.abs(mesh.divergence(psi)).max()*min(mesh.dr, mesh.dz)/bscale < 2e-14
    br, bz = mesh.cell_fields(.5*.03*r*r)
    assert_allclose(br, 0, atol=1e-15)
    assert_allclose(bz, .03, rtol=1e-14)


def test_stationary_uniform_state_axis_balance_and_conservation():
    mesh = MHDMesh(1, 1, 16, 24)
    r, _ = mesh.corners
    initial = from_primitive(mesh, .5*np.sqrt(MU0)*r*r, 1., 1.)
    solver = AxisymmetricMHD(mesh, initial)
    solver.advance(.1)
    assert_allclose(solver.state.conserved, initial.conserved, atol=3e-14, rtol=1e-13)
    assert_allclose(solver.state.psi, initial.psi, atol=1e-17)
    assert solver.diagnostics()["max_velocity_m_s"] < 1e-13


def test_full_one_dimensional_uct_dissipation():
    mesh = MHDMesh(1, 1, 12, 8)
    q = np.zeros((mesh.nr, mesh.nz, 6))
    q[..., 0] = 1; q[..., 3] = 1
    q[..., 1] = .1*np.sin(mesh.r[:, None]*5)
    q[..., 5] = np.sqrt(MU0)*(1+.2*np.cos(mesh.r[:, None]*3))
    br = np.zeros((mesh.nr+1, mesh.nz))
    bz = np.repeat(q[:, :1, 5], mesh.nz+1, axis=1)
    actual = ideal_emf(mesh, ghost_cells(q), br, bz, 5/3)
    left, right = q[:-1], q[1:]
    speed = np.maximum(np.abs(left[..., 1])+fast_speed(left, 0, 5/3), np.abs(right[..., 1])+fast_speed(right, 0, 5/3))
    expected = .5*(left[..., 1]*left[..., 5]+right[..., 1]*right[..., 5])-.5*speed*(right[..., 5]-left[..., 5])
    assert_allclose(actual[1:-1, :-1], expected, rtol=1e-14, atol=1e-18)


def test_cylindrical_fast_wave_converges_with_motion_and_advection():
    errors = []
    for n in (24, 48, 96):
        mesh = MHDMesh(1, 1, n, 4)
        initial, omega = radial_fast_mode(mesh)
        solver = AxisymmetricMHD(mesh, initial)
        t = np.pi/(2*omega)
        solver.advance(t)
        expected, _ = radial_fast_mode(mesh, time=t)
        q, exact = primitive(mesh, solver.state, 5/3), primitive(mesh, expected, 5/3)
        weights = mesh.volumes
        error = np.sqrt(np.sum((q[..., 1]-exact[..., 1])**2*weights)/np.sum(exact[..., 1]**2*weights))
        errors.append(error)
        assert q[..., 1].max() > 0
        assert abs(solver.diagnostics()["total_energy_relative_change"]) < 1e-12
        assert abs(solver.diagnostics()["mass_relative_change"]) < 1e-12
    assert errors[1] < .75*errors[0]
    assert errors[2] < .75*errors[1]
    assert errors[-1] < .08


def test_axial_advection_of_two_dimensional_flux_converges():
    errors = []
    for n in (12, 24, 48):
        mesh = MHDMesh(1, 1, n, 2*n)
        initial = advected_flux(mesh)
        solver = AxisymmetricMHD(mesh, initial)
        solver.advance(.15)
        exact = advected_flux(mesh, time=.15)
        errors.append(np.linalg.norm(solver.state.psi-exact.psi)/np.linalg.norm(exact.psi))
        field = max(np.abs(x).max() for x in mesh.face_fields(solver.state.psi))
        assert np.abs(mesh.divergence(solver.state.psi)).max()*mesh.dr/field < 5e-14
        assert abs(solver.diagnostics()["total_energy_relative_change"]) < 1e-12
    assert errors[1] < .75*errors[0]
    assert errors[2] < .75*errors[1]


def test_resistive_heating_and_decay_without_total_energy_injection():
    mesh = MHDMesh(1, 1, 64, 4)
    initial, k = resistive_bessel_mode(mesh)
    diffusivity, t = .5, .01
    solver = AxisymmetricMHD(mesh, initial, MHDParameters(resistivity=MU0*diffusivity))
    before = solver.diagnostics()
    solver.advance(t)
    after = solver.diagnostics()
    expected = initial.psi*np.exp(-diffusivity*k*k*t)
    assert np.linalg.norm(solver.state.psi-expected)/np.linalg.norm(expected) < .005
    magnetic_loss = before["magnetic_energy_j"]-after["magnetic_energy_j"]
    thermal_gain = after["internal_energy_j"]-before["internal_energy_j"]
    assert thermal_gain > 0 and after["physical_ohmic_heat_j"] > 0
    assert abs(after["total_energy_relative_change"]) < 1e-12
    assert_allclose(thermal_gain+after["kinetic_energy_j"], magnetic_loss, rtol=2e-7)
    assert_allclose(after["physical_ohmic_heat_j"], magnetic_loss, rtol=.03)


def test_sod_shock_tube_against_exact_solution():
    errors = []
    for n in (80, 160):
        mesh = MHDMesh(1, 1, 4, n)
        solver = AxisymmetricMHD(mesh, sod_initial(mesh), MHDParameters(gamma=1.4))
        solver.advance(.12)
        q = primitive(mesh, solver.state, 1.4)
        exact, _, _ = sod_exact(mesh.z, .12)
        central = abs(mesh.z) < .4
        errors.append(np.mean(abs(q[1, central, 0]-exact[central])))
        assert q[..., 3].min() > 0 and q[..., 0].min() > 0
        assert abs(solver.diagnostics()["total_energy_relative_change"]) < 1e-12
    assert errors[1] < .8*errors[0]
    assert errors[1] < .035


def test_nonphysical_states_and_unsafe_steps_are_rejected():
    mesh = MHDMesh(1, 1, 8, 8)
    flux = np.zeros((9, 9))
    with pytest.raises(NonphysicalState, match="pressure"):
        from_primitive(mesh, flux, 1, -1)
    with pytest.raises(NonphysicalState, match="density"):
        from_primitive(mesh, flux, -1, 1)
    state = from_primitive(mesh, flux, 1, 1)
    solver = AxisymmetricMHD(mesh, state)
    with pytest.raises(ValueError, match="Time step"):
        solver.step(2*solver.timestep_limit())
    assert_allclose(solver.state.conserved, state.conserved, rtol=0, atol=0)
    bad_flux = flux.copy(); bad_flux[-1] = 1e-5*np.cos(mesh.zf)
    bad = from_primitive(mesh, bad_flux, 1, 1)
    with pytest.raises(ValueError, match="radial wall"):
        AxisymmetricMHD(mesh, bad)


def test_rk2_temporal_order_at_fixed_spatial_grid():
    mesh = MHDMesh(1, 1, 24, 4)
    initial, omega = radial_fast_mode(mesh)
    end_time = np.pi/(2*omega)
    baseline = AxisymmetricMHD(mesh, initial).timestep_limit()
    velocities = []
    for fraction in (.8, .4, .2, .025):
        solver = AxisymmetricMHD(mesh, initial)
        steps = int(np.ceil(end_time/(fraction*baseline)))
        dt = end_time/steps
        for _ in range(steps):
            solver.step(dt)
        velocities.append(primitive(mesh, solver.state, 5/3)[..., 1])
    errors = [np.linalg.norm(v-velocities[-1]) for v in velocities[:-1]]
    assert errors[0]/errors[1] > 3.5
    assert errors[1]/errors[2] > 3.5


def test_magnetic_force_against_analytic_current_cross_field():
    errors = []
    for n in (24, 48, 96):
        mesh = MHDMesh(1, 1, n, 4)
        initial, b0, alpha = radial_pressure_balance(mesh, balanced=False)
        rhs, _ = AxisymmetricMHD(mesh, initial).rhs(initial)
        expected = radial_average(mesh, lambda r: 2*alpha*b0*b0/MU0*r*(1-alpha*r*r))
        actual = rhs.conserved[:, 0, 1]
        # Use the physical cylindrical volume norm, including the axis cells;
        # exclude only the wall region where a boundary stress is prescribed.
        interior = mesh.r < .9*mesh.radius
        weights = mesh.volumes[:, 0]
        errors.append(np.sqrt(np.sum((actual[interior]-expected[interior])**2*weights[interior])/
                              np.sum(expected[interior]**2*weights[interior])))
        assert np.all(actual[interior] > 0)
    assert errors[0]/errors[1] > 3.2
    assert errors[1]/errors[2] > 3.2


def test_equilibrium_numerical_flow_decreases_with_refinement():
    flows = []
    for n in (24, 48, 96):
        mesh = MHDMesh(1, 1, n, 4)
        initial, _, _ = radial_pressure_balance(mesh)
        solver = AxisymmetricMHD(mesh, initial)
        rhs, _ = solver.rhs(initial)
        assert np.abs(rhs.conserved[..., 1]).max() < 2e-12
        solver.advance(.2)
        flows.append(solver.diagnostics()["max_velocity_m_s"])
        assert abs(solver.diagnostics()["total_energy_relative_change"]) < 1e-12
    assert flows[1] < .8*flows[0]
    assert flows[2] < .8*flows[1]
