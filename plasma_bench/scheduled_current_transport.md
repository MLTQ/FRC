# Schedule independent gradient after stress workspace

`ScheduledCurrentTransport` inherits all released-workspace geometry, face and
advection methods. Its pair moves only velocity-gradient construction until
after the complete stress_rate evaluation and stress-array release. Each
arithmetic expression, density/state argument, coefficient and returned term
remains unchanged. The calls read the same immutable state and mapping.

Removing the moved assignment leaves identical statement ASTs, but the full
function's call order changes: do not claim the earlier del-only AST identity.
Require complete original/new coarse and fine signed physical workflows with
every input/rate/face/solver/independent-audit quantity bitwise equal, followed
by full fresh replays. Measure actual own-process peaks; theoretical array sizes
do not establish a resource bound.

No physical mass, density, boundary, basis, quadrature, residual or acceptance
budget changes. No diffusion, filtering, atmosphere, floor or energy repair.
Original sources/defaults remain frozen. This opt-in scheduling change does
not establish time accuracy, vacuum coupling, formation or confinement.
