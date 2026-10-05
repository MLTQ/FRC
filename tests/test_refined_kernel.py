"""Strict accuracy refinement and unchanged physical equations."""
import numpy as np
import pytest
from scipy.constants import elementary_charge as E, mu_0 as MU0
from examples.regular_plasma_reference import initialize
from plasma_bench.regular_plasma_step import field_from_packed
from plasma_bench.refined_regular_plasma import RefinedRegularPlasma


def localized_case():
    ref, d, u, original, field = initialize(dict(kind='wave', elements=1, points=8, field_nr=6, field_nz=8))
    u[...] = 0.; u[..., 0] = ref.number_density
    s = np.log(np.array([ref.electron_pressure, ref.ion_pressure])/(E*ref.number_density))
    s -= (d.thermal.gammas-1)*np.log(ref.number_density/d.thermal.n0)
    u[..., 4:] = ref.number_density*s
    b = np.zeros(original.geometry.hessian.shape[0]); b[56] = 1e-4*np.sqrt(MU0*ref.density)*ref.sound_speed/d.radius
    field = field_from_packed(original.geometry, b, ref.bias)
    return ref, u, original, field


def test_localized_magnetic_case_meets_every_original_residual_after_refinement():
    ref, u, original, field = localized_case(); before = u.copy()
    with pytest.raises(RuntimeError):
        original.derivative(u, field)
    refined = RefinedRegularPlasma(original.mapping, ref.ion_mass)
    result = refined.derivative(u, field)
    r = result.diagnostics['induction_residuals']
    assert result.residual < 1e-10
    assert r['outer_refinements'] > 0 and r['inner_tolerance'] < 1e-10
    assert r['full_faraday_backward'] < 1e-10
    assert np.array_equal(u, before)


def test_successful_original_case_is_bitwise_unchanged():
    ref, d, u, original, field = initialize(dict(kind='wave', elements=2))
    refined = RefinedRegularPlasma(original.mapping, ref.ion_mass)
    a, b = original.derivative(u, field), refined.derivative(u, field)
    for key in ('fluid', 'magnetic', 'current', 'electric'):
        assert np.array_equal(getattr(a, key), getattr(b, key))
    assert a.residual == b.residual
    assert b.diagnostics['induction_residuals']['outer_refinements'] == 0


def test_refinement_exhaustion_rejects_atomically():
    ref, u, original, field = localized_case(); before = u.copy()
    refined = RefinedRegularPlasma(original.mapping, ref.ion_mass, max_iterations=1)
    with pytest.raises(RuntimeError, match='strict inner refinement'):
        refined.derivative(u, field)
    assert np.array_equal(u, before)
