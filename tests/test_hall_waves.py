"""Independent species, cylindrical geometry, conservation and wave accuracy checks."""
import numpy as np
import pytest
from scipy.constants import elementary_charge as E, electron_mass as ME, proton_mass as MP, mu_0 as MU0
from plasma_bench.hall_parameters import HallParameters
from plasma_bench.hall_cylinder import HallCylinderMesh, LinearHallCylinder
from plasma_bench.hall_modes import LinearHallModes
from examples.hall_wave_reference import CylindricalWave, frequency


def setup(nr=24, nz=9, helicity=1, eta=0., model='two_fluid', length=2., radial_mode=1, axial_mode=1):
    p = HallParameters(1e18, MP, .001, eta, model)
    mesh = HallCylinderMesh(1., length, nr, nz)
    solver = LinearHallCylinder(mesh, p)
    wave = CylindricalWave(mesh.radius, mesh.length, p.number_density, p.ion_mass, p.bias,
                           eta, helicity, model, radial_mode, axial_mode)
    values = [x.real.copy() for x in wave.potentials(mesh.r[:, None], mesh.z[None, :])]
    for x in values:
        x[[0, -1]] = 0.
    return solver, wave, solver.from_nodal(*values)


@pytest.mark.parametrize('helicity', [-1, 1])
@pytest.mark.parametrize('eta', [0., 3e-5])
def test_parallel_reference_satisfies_original_species_momenta(helicity, eta):
    n, mi, bias, k = 1e18, MP, .001, 250.
    omega = frequency(n, mi, bias, eta, k, 0., helicity)
    b = np.array([1., 1j*helicity, 0.])*1e-8
    j = helicity*k/MU0*b
    v = -k*bias/(MU0*n*(mi+ME)*omega)*b
    vi = v+ME/(mi+ME)*j/(E*n)
    ve = v-mi/(mi+ME)*j/(E*n)
    electric = 1j*omega/(helicity*k)*b
    background = np.array([0., 0., bias])
    for mass, charge, velocity, friction in ((mi, E, vi, -E*n*eta*j), (ME, -E, ve, E*n*eta*j)):
        inertial = -1j*omega*mass*n*velocity
        electric_force = charge*n*electric
        magnetic_force = charge*n*np.cross(velocity, background)
        force = electric_force+magnetic_force+friction
        scale = sum(np.linalg.norm(term) for term in (inertial, electric_force, magnetic_force, friction))
        assert np.linalg.norm(inertial-force) < 3e-12*scale
    assert omega.imag <= 0


def test_parallel_high_k_limits_and_pair_mass():
    n, mi, bias = 1e18, MP, .001
    for h, mass in ((1, ME), (-1, mi)):
        omega = frequency(n, mi, bias, 0., 1e9, 0., h)
        assert omega.real == pytest.approx(E*bias/mass, rel=3e-11)
    pair = HallParameters(n, ME, bias)
    assert pair.hall_coefficient == 0
    assert pair.current_length_squared > 0
    assert frequency(n, ME, bias, 0., 200., 0., -1) == frequency(n, ME, bias, 0., 200., 0., 1)


def test_continuum_cylindrical_curl_and_divergence():
    _, wave, _ = setup(helicity=-1, eta=1e-5)
    r, z, delta = .41, .37, 2e-5
    b = wave.fields(r, z, complex_output=True)
    dr = (wave.fields(r+delta, z, complex_output=True)-wave.fields(r-delta, z, complex_output=True))/(2*delta)
    dz = (wave.fields(r, z+delta, complex_output=True)-wave.fields(r, z-delta, complex_output=True))/(2*delta)
    curl = np.array([-dz[1], dz[0]-dr[2], dr[1]+b[1]/r])
    assert np.linalg.norm(curl-wave.helicity*wave.total_k*b) < 2e-8*np.linalg.norm(curl)
    assert abs(dr[0]+b[0]/r+dz[2]) < 2e-8*np.linalg.norm(curl)
    axis = wave.fields(0., z)
    assert axis[0] == axis[1] == 0


def test_physical_volume_energy_heat_and_angular_momentum():
    solver, wave, state = setup(nr=12, eta=3e-5)
    # Add an independent smooth mean toroidal flow; it must retain angular momentum.
    mean = 2.*solver.mesh.r*(1-solver.mesh.r**2)
    extra = np.broadcast_to(mean[:, None], (13, 9)).copy()
    state += solver.from_nodal(np.zeros_like(extra), np.zeros_like(extra), np.zeros_like(extra), extra)
    r, z, weights = solver.quadrature_points(12, 27)
    b = solver.sample(state, r, z)
    v = solver.sample(state, r, z, 'velocity')
    j = solver.sample(state, r, z, 'current')
    physical = [np.sum(weights*np.sum(b*b, axis=-1))/(2*MU0),
                np.sum(weights*np.sum(v*v, axis=-1))*solver.parameters.density/2,
                np.sum(weights*np.sum(j*j, axis=-1))*MU0*solver.parameters.current_length_squared/2]
    assert list(solver.energies(state).values()) == pytest.approx(physical, rel=2e-12)
    initial = state.copy()
    step = solver.step(state, .2/wave.omega.real)
    middle_j = solver.sample(step.midpoint, r, z, 'current')
    actual_heat = solver.parameters.resistivity*np.sum(weights*np.sum(middle_j**2, axis=-1))*.2/wave.omega.real
    assert step.heat == pytest.approx(actual_heat, rel=2e-12)
    assert abs(step.balance_residual) < 2e-14*sum(physical)
    assert step.equation_relative_residual < 1e-15
    assert solver.angular_momentum(step.state) == pytest.approx(solver.angular_momentum(state), rel=2e-13)
    np.testing.assert_array_equal(state, initial)


@pytest.mark.parametrize('helicity', [-1, 1])
def test_actual_cylindrical_mesh_converges(helicity):
    errors = []
    for nr in (12, 24, 48):
        solver, wave, state = setup(nr=nr, helicity=helicity, eta=2e-5)
        time = .8/wave.omega.real
        for _ in range(80):
            state = solver.step(state, time/80).state
        r, z, weights = solver.quadrature_points(6)
        sample = solver.sample(state, r, z)
        exact = wave.fields(r, z, time)
        errors.append(np.sqrt(np.sum(weights*np.sum((sample-exact)**2, axis=-1))/np.sum(weights*np.sum(exact**2, axis=-1))))
    assert errors[0]/errors[1] > 1.8
    assert errors[1]/errors[2] > 1.8
    assert errors[-1] < .014


@pytest.mark.parametrize('helicity', [-1, 1])
def test_parallel_midpoint_time_order_phase_and_damping(helicity):
    p = HallParameters(1e18, MP, .001, 3e-5)
    k = 200.
    solver = LinearHallModes([k], [0.], p)
    omega = frequency(p.number_density, p.ion_mass, p.bias, p.resistivity, k, 0., helicity)
    ratio = -k*p.alfven/omega
    state = np.array([[[1., helicity, ratio, helicity*ratio]]], complex)*1e-8
    time = 2/omega.real
    errors = []
    for nsteps in (20, 40, 80):
        result = state.copy()
        heat = 0.
        for _ in range(nsteps):
            step = solver.step(result, time/nsteps)
            result, heat = step.state, heat+step.heat
        exact = state*np.exp(-1j*omega*time)
        errors.append(np.linalg.norm(result-exact)/np.linalg.norm(exact))
        balance = sum(solver.energies(result).values())+heat-sum(solver.energies(state).values())
        assert abs(balance) < 3e-14*sum(solver.energies(state).values())
    assert 3.95 < errors[0]/errors[1] < 4.05
    assert 3.95 < errors[1]/errors[2] < 4.05


def test_hall_generates_toroidal_field_and_torque_from_poloidal_perturbation():
    solver, wave, state = setup()
    state[..., 1:] = 0.
    dt = .1/wave.omega.real
    out = solver.step(state, dt).state
    assert np.linalg.norm(out[..., 1]) > .005*np.linalg.norm(state)
    out = solver.step(out, dt).state
    assert np.linalg.norm(out[..., 3]) > 1e-5*np.linalg.norm(state)
    control, _, initial = setup(model='mhd')
    initial[..., 1:] = 0.
    final = control.step(initial, dt).state
    np.testing.assert_array_equal(final[..., (1, 3)], 0.)


def test_axis_boundary_validation_and_stationary_state():
    solver, _, state = setup()
    state[:, 1:] = 0.
    np.testing.assert_array_equal(solver.step(state, 1.).state, state)
    for bad in (0., -1., np.nan):
        with pytest.raises(ValueError):
            solver.step(state, bad)
    with pytest.raises(ValueError):
        HallCylinderMesh(1., 1., 8, 8)
    with pytest.raises(ValueError):
        HallParameters(0., MP, .01)
    broken = state.copy()
    broken[0, 1, 0] = 1.
    with pytest.raises(ValueError):
        solver.step(broken, .01)
    values = solver.sample(state, np.array([0., 1.]), 0.)
    np.testing.assert_allclose(values[:, :2], 0., atol=1e-25)


def test_large_midpoint_step_has_bad_phase_despite_exact_energy():
    p = HallParameters(1e18, MP, .001, model='hall')
    k, h = 20., 1
    omega = frequency(p.number_density, MP, p.bias, 0., k, 0., h, 'hall').real
    solver = LinearHallModes([k], [0.], p)
    state = np.array([[[1., h, -k*p.alfven/omega, -h*k*p.alfven/omega]]], complex)*1e-8
    step = solver.step(state, 2/omega)
    assert abs(step.balance_residual) < 1e-14*sum(solver.energies(state).values())
    phase = -np.angle(step.state[0, 0, 0]/state[0, 0, 0])
    assert phase == pytest.approx(np.pi/2, rel=2e-14)
    assert abs(phase/2-1) > .2


def test_exact_time_oracle_and_local_continuum_heat():
    from scipy.integrate import quad
    from examples.hall_time_reference import exact_spatial_evolution
    solver, wave, state = setup(nr=8, eta=3e-5)
    time = 2/wave.omega.real
    exact = exact_spatial_evolution(solver, state, time)
    errors = []
    for count in (40, 80, 160):
        result = state.copy()
        for _ in range(count):
            result = solver.step(result, time/count).state
        errors.append(np.linalg.norm(result-exact))
    assert errors[0]/errors[1] > 3.8
    assert errors[1]/errors[2] > 3.8
    r, z = .3, .43
    numerical_heat = quad(lambda t: wave.resistivity*np.sum(wave.fields(r, z, t, 'current')**2), 0, time, epsabs=1e-30)[0]
    assert wave.heat_density(r, z, time) == pytest.approx(numerical_heat, rel=1e-12)


def test_signed_bias_reversal_exchanges_helicities():
    n, bias, k = 1e18, .001, 100.
    for h in (-1, 1):
        positive = frequency(n, MP, bias, 3e-5, k, 0., h)
        negative = frequency(n, MP, -bias, 3e-5, k, 0., -h)
        assert positive == negative
        p = HallParameters(n, MP, -bias, 3e-5)
        solver = LinearHallModes([k], [0.], p)
        ratio = -k*p.alfven/negative
        state = np.array([[[1., -h, ratio, -h*ratio]]], complex)*1e-8
        step = solver.step(state, .01/negative.real)
        exact = state*np.exp(-1j*negative*.01/negative.real)
        assert np.linalg.norm(step.state-exact)/np.linalg.norm(exact) < 1e-7
