"""Full physical uniform-inertia inverse and variable-density equivalence."""
import numpy as np
from scipy.constants import elementary_charge as E
from examples.strong_wall_reference import initialize
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.wall_induction_core import WallInductionCore
from plasma_bench.uniform_wall_metric import UniformWallMetric
from plasma_bench.uniform_wall_plasma import preconditioned
from plasma_bench.modal_wall_coordinates import ModalWallCoordinates


def test_uniform_metric_solves_full_independent_physical_operator_and_is_spd():
    _, _, u, s, b = initialize(dict(kind='contact', elements=1, points=8, field_nr=6, field_nz=8))
    c = ModalWallCoordinates(s, u, b); u = c.background()
    model = RegularCurrent(s.mapping, u[..., 0], s.ion_mass); core = WallInductionCore(model)
    inverse = UniformWallMetric(model, core.lift)
    vectors = np.random.default_rng(884).normal(size=(3, core.lift.shape[1]))
    solved = []
    for value in vectors:
        x = inverse.solve(value); magnetic = core.lift@x; j = model.current(magnetic)
        dual = s.geometry.hessian@magnetic-s.geometry.hessian@s.geometry.faraday@core.mass_inverse(model.inertia@j)
        physical = core.lift.T@dual
        assert np.linalg.norm(physical-value)/np.linalg.norm(value) < 1e-10
        assert value@x > 0
        solved.append(x)
    assert abs(vectors[0]@solved[1]-vectors[1]@solved[0])/(np.linalg.norm(vectors[0])*np.linalg.norm(solved[1])) < 1e-13


def test_actual_nonuniform_rhs_remains_equivalent_under_strict_bounds():
    _, _, u, original, b = initialize(dict(kind='contact', elements=1, points=8, field_nr=6, field_nz=8))
    c = ModalWallCoordinates(original, u, b); changed = preconditioned(original)
    for scale in (1e-4, 2e-4):
        du, db = c.unpack(scale*np.random.default_rng(504).normal(size=c.size))
        state = u+du; field = c.field(db)
        a = original.derivative(state, field); result = changed.derivative(state, field)
        assert np.array_equal(result.fluid, a.fluid)
        for key in ('magnetic', 'current', 'electric'):
            x, y = getattr(result, key), getattr(a, key)
            assert np.linalg.norm(x-y)/np.linalg.norm(y) < 1e-10
        assert result.residual <= original.solver_options['tolerance']
