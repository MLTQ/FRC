"""Independent linearization, original collocation and atomic-failure controls."""
import numpy as np
import pytest
from examples.regular_plasma_reference import initialize
from examples.regular_plasma_linear import StationaryLinearPlasma
from plasma_bench.refined_regular_plasma import RefinedRegularPlasma
from plasma_bench.plasma_coordinates import PlasmaCoordinates
from plasma_bench.plasma_preconditioner import PlasmaPreconditioner
from plasma_bench.complex_plasma_inverse import ComplexPlasmaInverse
from plasma_bench.implicit_plasma_midpoint import ImplicitPlasmaMidpoint
from plasma_bench.implicit_plasma_gauss import ImplicitPlasmaGauss
from plasma_bench.implicit_plasma_radau import ImplicitPlasmaRadau
from plasma_bench.complex_radau_inverse import ComplexRadauInverse
from plasma_bench.regular_fluid_linear import operator


def setup(kind='wave', **extra):
    ref, d, u, old, b = initialize(dict(kind=kind, elements=1, points=8,
                                      field_nr=6, field_nz=8, **extra))
    return ref, d, u, RefinedRegularPlasma(old.mapping, ref.ion_mass, **old.solver_options), b


def test_uniform_moving_fluid_linearization_converges():
    _, d, u, system, b = setup('contact')
    c = PlasmaCoordinates(system, u, b); u = c.background()
    direction, _ = c.unpack(np.random.default_rng(875).normal(size=c.size))
    expected = (operator(d, u)@direction.ravel()).reshape(d.shape)
    errors = []
    for epsilon in (1e-3, 1e-4):
        difference = (d.derivative(u+epsilon*direction)[0]-d.derivative(u-epsilon*direction)[0])/(2*epsilon)
        scaled = (difference-expected).ravel()[c.fluid_indices]/c.fluid_scale
        errors.append(np.linalg.norm(scaled)/np.linalg.norm(expected.ravel()[c.fluid_indices]/c.fluid_scale))
    assert 95 < errors[0]/errors[1] < 105
    assert errors[1] < 1e-6


@pytest.mark.parametrize('method', ['gauss', 'radau'])
def test_physical_preconditioner_matches_independent_rest_and_stage_equations(method):
    _, _, u, system, b = setup()
    u[..., 1:4] = 0.
    c = PlasmaCoordinates(system, u, b); background = c.background()
    p = PlasmaPreconditioner(c); reference = StationaryLinearPlasma(system, background, b.axial_bias)
    matrix = np.block([[p.fluid.toarray(), p.to_fluid], [p.to_magnetic, p.magnetic]])
    for direction in np.random.default_rng(91).normal(size=(3, c.size)):
        du, db = c.unpack(direction)
        fluid, magnetic, _ = reference.derivative(du, db)
        expected = c.pack(fluid, magnetic)
        assert np.linalg.norm(matrix@direction-expected)/np.linalg.norm(expected) < system.solver_options['tolerance']
    dt = 1e-6
    if method == 'gauss':
        inverse = ComplexPlasmaInverse(p); inverse.factor(dt*(.25+1j*np.sqrt(3)/12))
        a = np.array([[.25, .25-np.sqrt(3)/6], [.25+np.sqrt(3)/6, .25]])
    else:
        inverse = ComplexRadauInverse(p); inverse.factor(dt*(1/3+1j*np.sqrt(2)/6))
        a = np.array([[5/12, -1/12], [3/4, 1/4]])
    residual = np.random.default_rng(85).normal(size=(2, c.size))
    correction = inverse.stages(residual)
    error = correction-dt*a@(correction@matrix.T)-residual
    assert np.linalg.norm(error)/np.linalg.norm(residual) < 1e-10


@pytest.mark.parametrize('method', [ImplicitPlasmaMidpoint, ImplicitPlasmaGauss, ImplicitPlasmaRadau])
def test_original_time_equations_and_atomic_exhaustion(method):
    ref, d, u, system, b = setup()
    before = u.copy(); magnetic = b.packed_increment(system.geometry).copy()
    solver = method(system, u, b); c = solver.coordinates; dt = 2*np.pi/abs(ref.omega)/120
    evolved, field, diagnostic = solver.step(u, b, dt)
    change = c.pack(evolved-u, field.packed_increment(system.geometry)-magnetic)
    if method is ImplicitPlasmaMidpoint:
        rate = system.derivative(.5*(u+evolved), c.field(.5*(magnetic+field.packed_increment(system.geometry))))
        error = change-dt*c.pack(rate.fluid, rate.magnetic)
    else:
        rates = []; actual = []
        for y in diagnostic['stage_changes']:
            du, db = c.unpack(y); state = u+du; packed = magnetic+db
            actual.append(c.pack(state-u, packed-magnetic))
            rate = system.derivative(state, c.field(packed)); rates.append(c.pack(rate.fluid, rate.magnetic))
        errors = np.asarray(actual)-dt*solver.a@rates
        weights = solver.weights if method is ImplicitPlasmaRadau else np.array([.5, .5])
        error = np.r_[errors.ravel(), change-dt*weights@rates]
        if method is ImplicitPlasmaRadau:
            du, _ = c.unpack(diagnostic['stage_changes'][1])
            assert np.array_equal(evolved, u+du)
    assert np.max(abs(error)) <= diagnostic['residual_limit']
    assert np.all(evolved[..., 0] > 0)
    assert np.array_equal(u, before) and np.array_equal(b.packed_increment(system.geometry), magnetic)
    ref, _, u, system, b = setup('contact')
    before = u.copy(); magnetic = b.packed_increment(system.geometry).copy()
    solver = method(system, u, b, max_iterations=1)
    with pytest.raises(RuntimeError, match='did not converge'):
        solver.step(u, b, ref.period/120)
    with pytest.raises(ValueError): solver.step(u, b, -dt)
    assert np.array_equal(u, before) and np.array_equal(b.packed_increment(system.geometry), magnetic)


def test_saved_contact_chord_failure_recovers_with_unchanged_bound():
    from pathlib import Path
    ref, _, u, old, b = initialize(dict(kind='contact', elements=4))
    system = RefinedRegularPlasma(old.mapping, ref.ion_mass, **old.solver_options)
    solver = ImplicitPlasmaGauss(system, u, b)
    path = Path(__file__).resolve().parents[1]/'outputs/implicit_gauss_contact_pilot/progress.npz'
    with np.load(path) as raw:
        before = raw['state'].copy(); magnetic = raw['magnetic'].copy()
    state = before.copy(); field = solver.coordinates.field(magnetic)
    _, _, diagnostic = solver.step(state, field, ref.period/120)
    assert diagnostic['newton_corrections'] >= 1
    assert diagnostic['residual'] < 1e-12
    assert diagnostic['residual'] <= diagnostic['residual_limit']
    assert diagnostic['max_induction_residual'] <= 1e-10
    assert np.array_equal(state, before)
    assert np.array_equal(field.packed_increment(system.geometry), magnetic)
