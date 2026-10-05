"""Independent full-parent checks for the bounded moving plasma assembly."""
import numpy as np
import pytest
from examples.regular_plasma_reference import initialize
from examples.regular_plasma_rates import run_case
from examples.regular_species_reference import initialize as heterogeneous
from plasma_bench.regular_plasma import RegularPlasma
from plasma_bench.regular_plasma_step import rk4


@pytest.mark.parametrize('kind', ['wave', 'contact'])
def test_full_original_parent_rates_refine(kind):
    a, _ = run_case(dict(kind=kind, elements=4))
    b, _ = run_case(dict(kind=kind, elements=8))
    for key in ('density_rate', 'electron_pressure_rate', 'ion_pressure_rate', 'electric'):
        coarse = a['independent']['norms'][key]['relative_L2']
        fine = b['independent']['norms'][key]['relative_L2']
        assert coarse/fine > 15, (kind, key, coarse, fine)
        assert fine < 6e-4
    if kind == 'wave':
        for key in ('velocity_rate', 'magnetic_rate', 'current_rate'):
            assert b['independent']['norms'][key]['relative_L2'] < 1e-4
    assert b['residual'] < 1e-10


def test_full_heterogeneous_energy_and_species_entropy():
    ref, d, u, model, field = heterogeneous(dict(elements=4))
    system = RegularPlasma(model.mapping, ref.ion_mass)
    result = system.derivative(u, field)
    powers = result.diagnostics['powers']
    assert abs(sum(powers.values()))/sum(abs(p) for p in powers.values()) < 1e-10
    for k in (0, 4, 5):
        integrated = d.volume*result.fluid[..., k]
        assert abs(integrated.sum())/abs(integrated).sum() < 1e-12


def test_geometry_cache_is_bitwise_with_changing_density():
    ref, d, u, system, field = initialize(dict(kind='wave', elements=2))
    system.derivative(u, field)
    changed = u.copy()
    changed[..., 0] *= 1+.02*(d.r/d.radius)**2
    changed[..., 4:6] += .03*changed[..., 0, None]
    cached = system.derivative(changed, field)
    fresh = RegularPlasma(system.mapping, ref.ion_mass, cache_geometry=False).derivative(changed, field)
    for key in ('fluid', 'magnetic', 'current', 'electric'):
        assert np.array_equal(getattr(cached, key), getattr(fresh, key))
    assert cached.diagnostics == fresh.diagnostics


def test_coupled_failures_are_atomic():
    ref, d, u, system, field = initialize(dict(kind='wave', elements=2))
    before = u.copy(); magnetic = field.packed_increment(system.geometry).copy()
    with pytest.raises(ValueError):
        rk4(system, u, field, -1.)
    bad = u.copy(); bad[1, 1, 0, 0] = -1.
    with pytest.raises(ValueError):
        system.derivative(bad, field)
    with pytest.raises(RuntimeError):
        RegularPlasma(system.mapping, ref.ion_mass, max_iterations=1).derivative(u, field)
    assert np.array_equal(u, before)
    assert np.array_equal(field.packed_increment(system.geometry), magnetic)


def test_full_explicit_step_preserves_inputs_and_axis_contract():
    ref, d, u, system, field = initialize(dict(kind='wave', elements=2))
    before = u.copy(); magnetic = field.packed_increment(system.geometry).copy()
    evolved, b, diagnostics = rk4(system, u, field, 2*np.pi/ref.omega/120)
    assert np.array_equal(u, before)
    assert np.array_equal(field.packed_increment(system.geometry), magnetic)
    assert np.count_nonzero(evolved[0, 0, :, 1:3]) == 0
    assert np.all(evolved[..., 0] > 0)
    assert diagnostics['max_residual'] < 1e-10
    assert not np.array_equal(b.packed_increment(system.geometry), magnetic)
