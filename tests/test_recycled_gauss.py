"""Independent actual-stage checks and transactional continuation memory."""
import copy
import numpy as np
import pytest
from examples.strong_wall_reference import initialize
from plasma_bench.continued_modal_gauss import ContinuedModalGauss
from plasma_bench.recycled_modal_gauss import RecycledModalGauss
from plasma_bench.guarded_modal_gauss import GuardedModalGauss


@pytest.mark.parametrize('solver_class', [ContinuedModalGauss, RecycledModalGauss, GuardedModalGauss])
def test_continued_actual_stages_and_failed_step_preserve_previous_history(solver_class):
    ref, _, state, system, field = initialize(dict(kind='contact', elements=1, points=8, field_nr=6, field_nz=8))
    solver = solver_class(system, state, field); c = solver.coordinates; dt = ref.period/120
    for index in range(2):
        before = state.copy(); packed = field.packed_increment(system.geometry).copy()
        updated, evolved, record = solver.step(state, field, dt)
        assert record['continued'] == bool(index)
        actual = []; rates = []
        for change in record['stage_changes']:
            du, db = c.unpack(change); stage = state+du; magnetic = packed+db
            rate = system.derivative(stage, c.field(magnetic))
            rates.append(c.pack(rate.fluid, rate.magnetic))
            actual.append(c.pack(stage-state, magnetic-packed))
        residual = np.array(actual)-dt*solver.a@rates
        endpoint = c.pack(updated-state, evolved.packed_increment(system.geometry)-packed)-.5*dt*np.sum(rates, axis=0)
        assert max(np.max(abs(residual)), np.max(abs(endpoint))) <= record['residual_limit']
        assert np.array_equal(state, before)
        assert np.array_equal(field.packed_increment(system.geometry), packed)
        state, field = updated, evolved
    memory = copy.deepcopy(solver._continuation); before = state.copy()
    packed = field.packed_increment(system.geometry).copy(); solver.max_iterations = 1
    with pytest.raises(RuntimeError): solver.step(state, field, dt*17)
    assert np.array_equal(state, before)
    assert np.array_equal(field.packed_increment(system.geometry), packed)
    for key, value in memory.items():
        assert np.array_equal(solver._continuation[key], value)
