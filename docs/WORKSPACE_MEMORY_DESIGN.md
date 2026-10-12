# Reduce temporary storage before independent field refinement

Continue from5a60b36, classified as verified progress. All preceding workers
are terminal. The full FRC goal stays open; no unit tests are added or run.

The first candidate changes only local object lifetimes in the original
current transport and electromagnetic functions. Release quadrature gradient,
stress, velocity, current and field temporaries after their last use. Preserve
every arithmetic expression, call order, coefficient, return term, boundary,
physical mass/density, residual and acceptance budget. The original functions
and all archived sources remain untouched. New subclasses are opt-in.

Require AST identity after removing only added del statements, then complete
coarse and fine signed physical workflows against the verified tensor baseline.
Every input, derivative, physical/face fingerprint and non-runtime metric must
agree bitwise; also perform fresh full candidate replays. This source check is
supporting evidence, not a substitute for the independent physical workflows.
Record own-process peak RSS and disk headroom; actual reductions are measured,
not inferred by subtracting theoretical array sizes or overlapping intervals.

The lifetime-only full fine calculation matches every original coefficient,
solver metric and independently audited physical quantity bitwise, but peaks
at6.873GiB. Its fresh full replay is bitwise identical but peaks7.430GiB;
these observations do not establish reliable memory reduction. A second candidate retains only
geometry data in cached templates after returning fully rebound operator copies.
It collects unreachable model/inertia reference cycles after the full parent
derivative returns. Actual current-state copies remain attached during work;
all original approximate metric values and physical operators are preserved.
This additional lifetime change needs its own full coarse/fine comparisons,
resource measurement and complete fresh replays. Coarse evidence cannot prove
the final fine peak or doubled-grid feasibility.

Only attempt larger axial/radial field grids when measured headroom justifies
them. Keep the original8GiB worker/256MiB disk resource guard with bounded
self-termination. Resource interruption is distinct from numerical rejection.
If scratch lifetime changes are insufficient, stream exact quadrature work in
a separately verified variant; never alter physical equations to fit resources.

The detached-cache fine calculation reproduces the released-workspace baseline
bitwise and peaks6.446GiB; its full replay completes at7.056GiB with identical
scientific results. Proceed with axial field cells256→288 at fixed fluid64/256 and
radial field256. This12.5% increase is an explicit resource-bounded refinement
step, not a substitute for the pending doubled grids or a convergence-order
claim. Simple storage scaling is planning evidence only; measure the real
guarded run. Require both signed solves, independent complete volume/face and
current-component/strong-curl comparisons, then a full fresh replay.

This step is now complete: both signed rates and every physical/face array
replay bitwise, with independent q16/20 component/strong-curl audits also fully
replayed. Jr error falls5.9694%→3.0532%, strong-curl Jr error5.9425%→3.0338%.
The weak/strong difference slightly increases, so the improvement concerns
the actual magnetic response, not a repaired weak Ampere identity. Axial current
error remains1.1140%. See WORKSPACE_MEMORY_VERIFICATION.md for scope, resources,
work defects and unfinished refinement/amplitude/time/formation requirements.

All mathematical assumptions remain the original filled periodic cylinder,
impermeable radial fluid wall, prescribed constant bias and bounded compatible
field space. There is no added diffusion, filtering, atmosphere, artificial
mass, density floor or energy repair. These instantaneous fixtures cannot prove
formation, time convergence, exact moving vacuum, 3-D stability or confinement.
