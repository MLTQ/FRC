"""Independent original-rate, invariant, trace and atomic-failure checks."""
import numpy as np
import pytest
from examples.regular_species_reference import initialize
from examples.regular_species_study import run_case
from plasma_bench.regular_species import RegularSpecies
from plasma_bench.regular_species_transport import FrozenSpeciesTransport
from plasma_bench.regular_current import RegularCurrent


def test_original_species_rates_refine_with_independent_temperatures():
    coarse, _ = run_case(dict(elements=4))
    fine, _ = run_case(dict(elements=8))
    for k in range(2):
        assert coarse['species'][k]['point_rate']['relative_L2']/fine['species'][k]['point_rate']['relative_L2'] > 10
        assert fine['species'][k]['point_rate']['relative_L2'] < .002
        assert fine['species'][k]['entropy_backward'] < 1e-12
        assert fine['species'][k]['quadratic_backward'] < 1e-12
        assert fine['species'][k]['work_backward'] < 1e-12
    assert fine['total_grid_power_backward'] < 1e-10
    assert fine['controls']['point_skew'] > 1


def test_discontinuous_traces_conserve_entropy_and_pair_actual_grid_work():
    _, d, u, m, field = initialize(dict(elements=3, field_nr=12))
    rng = np.random.default_rng(2937)
    u[..., 4:6] += u[..., 0, None]*rng.uniform(-.2, .2, u[..., 4:6].shape)
    before = u.copy()
    pair = RegularSpecies(m).pair(u, field.current(m))
    f = d.fields(u)
    for k in range(2):
        rate = pair['fluid_rate'][..., k+4]
        assert abs(np.sum(d.volume*rate))/np.sum(d.volume*abs(rate)) < 1e-12
        q = d.volume*f['entropy'][..., k]*rate
        assert abs(q.sum())/abs(q).sum() < 1e-12
        h = f['specific_internal'][..., k]
        work = np.sum(d.volume*h*rate)
        current_work = field.current(m)@pair['species_forcing'][k]
        assert abs(work-current_work) < 1e-11*np.sum(d.volume*abs(h*rate))
    assert np.count_nonzero(pair['fluid_rate'][..., :4]) == 0
    assert np.array_equal(u, before)


def test_frozen_operator_is_skew_in_actual_number_weighted_metric():
    _, d, u, m, field = initialize(dict(elements=1, points=8, field_nz=8))
    transport = FrozenSpeciesTransport(m, field.current(m))
    root = np.sqrt(d.volume*u[..., 0])
    size = root.size
    matrix = np.column_stack([(root*transport.rate(np.eye(size)[i].reshape(root.shape)/root)).ravel()
                              for i in range(size)])
    assert np.linalg.norm(matrix+matrix.T)/np.linalg.norm(matrix) < 2e-12
    eigenvalues = np.linalg.eigvals(matrix)
    assert np.max(abs(eigenvalues.real))/np.max(abs(eigenvalues)) < 2e-12
    assert np.linalg.norm(transport.rate(np.ones(root.shape)))/np.linalg.norm(matrix) < 2e-12


def test_species_mismatch_overflow_and_bad_step_reject_without_mutation():
    _, d, u, m, field = initialize(dict(elements=2))
    before = u.copy()
    model = RegularSpecies(m)
    bad = u.copy()
    bad[..., 0] *= 1.01
    with pytest.raises(ValueError):
        model.pair(bad, field.current(m))
    bad = u.copy()
    bad[..., 4] = 1000*bad[..., 0]
    with pytest.raises(ValueError):
        model.pair(bad, field.current(m))
    transport = FrozenSpeciesTransport(m, field.current(m))
    s = u[..., 4]/u[..., 0]
    with pytest.raises(ValueError):
        transport.step(s, -1)
    assert np.array_equal(u, before)


def test_uniform_temperatures_remain_in_analytic_pressure_curl_kernel():
    _, d, u, m, field = initialize(dict(elements=4))
    for k, gamma in enumerate(d.thermal.gammas):
        u[..., k+4] = u[..., 0]*(np.log(2.+k)-(gamma-1)*np.log(u[..., 0]/d.thermal.n0))
    pair = RegularSpecies(m).pair(u, field.current(m))
    assert np.linalg.norm(pair['reduced_forcing'])/np.linalg.norm(pair['forcing']) < 2e-12


def test_numeric_input_dtypes_do_not_truncate_species_rates():
    ref, d, u, m, field = initialize(dict(elements=2))
    integer = np.zeros(d.shape, dtype=np.int64)
    integer[..., 0] = 10**18
    integer[..., 4] = np.rint(u[..., 4]/u[..., 0]*10**18).astype(np.int64)
    integer[..., 5] = np.rint(u[..., 5]/u[..., 0]*10**18).astype(np.int64)
    model = RegularCurrent(m.mapping, integer[..., 0], ref.ion_mass)
    pair = RegularSpecies(model)
    current = field.current(model)
    result = pair.pair(integer, current)['fluid_rate']
    assert result.dtype == np.float64
    assert np.array_equal(result, pair.pair(integer.astype(float), current)['fluid_rate'])
