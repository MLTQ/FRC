"""Observation-only iteration must preserve original accepted numerical values."""
import numpy as np
from examples.strong_wall_reference import initialize
from examples.wall_roundoff_core import WallRoundoffCore
from plasma_bench.balanced_wall_core import BalancedWallCore
from plasma_bench.uniform_wall_plasma import UniformWallInduction
from plasma_bench.regular_current import RegularCurrent


def test_residual_observations_leave_original_iteration_bitwise_unchanged():
    _,_,u,s,_=initialize(dict(kind='contact',elements=2,points=8,field_nr=6,field_nz=8))
    model=RegularCurrent(s.mapping,u[...,0],s.ion_mass);solver=UniformWallInduction(model,tolerance=1e-10)
    force=np.random.default_rng(389).normal(size=model.mass.shape[0]);before=force.copy()
    original=BalancedWallCore.solve(solver,force);observed=WallRoundoffCore.solve(solver,force)
    for key in ('magnetic','current','electric'):assert np.array_equal(getattr(original,key),getattr(observed,key))
    assert original.residual_components==observed.residual_components
    assert original.iterations==observed.iterations and original.refinements==observed.refinements
    assert np.array_equal(force,before)
    assert all(np.isfinite(v).all() for v in solver.roundoff_arrays.values())
