"""Uniform physical symbols, full Gauss equations and modal atomic rejection."""
import numpy as np
import pytest
from examples.strong_wall_reference import initialize
from plasma_bench.strong_wall_coordinates import StrongWallCoordinates
from plasma_bench.strong_wall_preconditioner import StrongWallPreconditioner
from plasma_bench.modal_wall_coordinates import ModalWallCoordinates
from plasma_bench.modal_wall_inverse import ModalWallInverse
from plasma_bench.modal_wall_gauss import ModalWallGauss
from plasma_bench.accelerated_modal_gauss import AcceleratedModalGauss


def action(inverse, value):
    c = inverse.coordinates; modes = c.fourier(value); result = np.empty_like(modes)
    n = c.fluid_row_count
    for k, (f, fb, bf, b) in enumerate(inverse.operators):
        result[:, k] = np.r_[f@modes[:n, k]+fb@modes[n:, k], bf@modes[:n, k]+b@modes[n:, k]]
    return c.physical(result)


@pytest.mark.parametrize('moving', [False, True])
def test_modal_symbols_match_full_physical_blocks_and_complex_stage_equations(moving):
    _, _, u, system, field = initialize(dict(kind='contact', elements=1, points=8, field_nr=6, field_nz=8))
    if not moving: u[..., 1:4] = 0.
    c = ModalWallCoordinates(system, u, field); inverse = ModalWallInverse(c)
    old = StrongWallCoordinates(system, u, field); p = StrongWallPreconditioner(old)
    matrix = np.block([[p.fluid.toarray(), p.to_fluid], [p.to_magnetic, p.magnetic]])
    # Random fields contain every Fourier mode, including zero and Nyquist.
    for direction in np.random.default_rng(118).normal(size=(3, c.size)):
        physical = c.unpack(direction); expected = c.pack(*old.unpack(matrix@old.pack(*physical)))
        value = action(inverse, direction)
        assert np.linalg.norm(value-expected)/np.linalg.norm(expected) < 1e-10
        assert np.linalg.norm(value.imag)/np.linalg.norm(value) < 1e-12
        assert np.max(abs(c.pack(*c.unpack(direction))-direction)) < 2e-14
    dt = 1e-6; inverse.factor(dt*(.25+1j*np.sqrt(3)/12))
    rhs = np.random.default_rng(97).normal(size=(2, c.size))
    correction = inverse.stages(rhs)
    a = np.array([[.25, .25-np.sqrt(3)/6], [.25+np.sqrt(3)/6, .25]])
    residual = correction-dt*a@np.array([action(inverse, v) for v in correction])-rhs
    assert np.linalg.norm(residual)/np.linalg.norm(rhs) < 1e-10


def test_uniform_nonlinear_direction_converges_to_modal_parent_operator():
    _, _, u, system, field = initialize(dict(kind='contact', elements=1, points=8, field_nr=6, field_nz=8))
    u[..., 1:4] = 0.; c = ModalWallCoordinates(system, u, field); u = c.background()
    inverse = ModalWallInverse(c); direction = np.random.default_rng(61).normal(size=c.size)
    du, db = c.unpack(direction); expected = action(inverse, direction).real; errors = []
    for epsilon in (1e-5, 1e-6):
        positive = system.derivative(u+epsilon*du, c.field(epsilon*db))
        negative = system.derivative(u-epsilon*du, c.field(-epsilon*db))
        value = c.pack(positive.fluid-negative.fluid, positive.magnetic-negative.magnetic)/(2*epsilon)
        errors.append(np.linalg.norm(value-expected)/np.linalg.norm(expected))
    assert errors[1] < 1e-8
    assert 80 < errors[0]/errors[1] < 120


@pytest.mark.parametrize("solver_class", [ModalWallGauss, AcceleratedModalGauss])
def test_actual_modal_gauss_stages_and_atomic_failure(solver_class):
    ref, _, u, s, b = initialize(dict(kind='wave', elements=1, points=8, field_nr=6, field_nz=8))
    solver = solver_class(s, u, b); c = solver.coordinates; dt = 2*np.pi/abs(ref.omega)/60
    before = u.copy(); packed = b.packed_increment(s.geometry)
    evolved, field, diagnostic = solver.step(u, b, dt)
    rates = []; stages = []
    for change in diagnostic['stage_changes']:
        du, db = c.unpack(change); stage = u+du; magnetic = packed+db
        rate = s.derivative(stage, c.field(magnetic)); rates.append(c.pack(rate.fluid, rate.magnetic))
        stages.append(c.pack(stage-u, magnetic-packed))
        assert np.all(stage[-1, -1, :, 1] == 0)
    residual = np.array(stages)-dt*solver.a@rates
    endpoint = c.pack(evolved-u, field.packed_increment(s.geometry)-packed)-.5*dt*np.sum(rates, axis=0)
    assert max(np.max(abs(residual)), np.max(abs(endpoint))) <= diagnostic['residual_limit']
    assert np.array_equal(u, before) and np.array_equal(b.packed_increment(s.geometry), packed)
    ref, _, u, s, b = initialize(dict(kind='contact', elements=1, points=8, field_nr=6, field_nz=8))
    before = u.copy(); packed = b.packed_increment(s.geometry)
    exhausted = solver_class(s, u, b, max_iterations=1)
    with pytest.raises(RuntimeError, match='did not converge'): exhausted.step(u, b, ref.period/120)
    assert np.array_equal(u, before) and np.array_equal(b.packed_increment(s.geometry), packed)
