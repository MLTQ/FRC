# Geometry templates without retained completed density models

`DetachedMetricWallPlasma.operators` first obtains the original fully rebound
operator copies. Only then it clears model references in the cached templates.
The copies used for the actual evaluation continue to hold the complete current
density/mass/inertia model. All differential maps, physical matrices, lift,
species projection and approximate metric/factor values remain unchanged.

This relies on the inspected original cache contract: CurrentTransport retains
geometry/differential/face maps; RegularSpecies uses a SpeciesProjection holding
the geometry/fluid mapping; the induction approximate inverse stores numeric
factors and reference density. None needs a template's old model after rebinding.
Returned copies are never detached while computing a rate. With caching disabled,
the parent path is unchanged.

`derivative` evaluates the complete parent method including strict checks and
work diagnostics, then invokes garbage collection after its frame has returned.
This releases unreachable model/inertia cycles formed by LinearOperator bound
methods. The returned rates/diagnostics and caller inputs remain referenced.
Collection is workspace management, not a physical operation or energy repair.

`detached_metric` constructs this opt-in system with the exact original settings.
Require complete coarse/fine signed old/new physical comparisons and fresh
replays. Actual output must be bitwise identical; source reasoning alone does
not certify equivalence. Record own-worker memory and elapsed time. No unit tests,
physical-equation changes or default time-integrator adoption are introduced.
