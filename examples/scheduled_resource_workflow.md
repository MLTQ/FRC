# Stage memory observations during a full scheduled replay

`run(destination, observed)` repeats the entire scheduled F64/256, B256/256
signed calculation, including q16/20 volume and face audits. It requires every
saved coefficient and non-runtime scientific metric to reproduce the unobserved
`scheduled_current_fine` baseline. No unit tests are run.

Entry/exit observations wrap the actually dispatched scheduled transport,
released advection/electromagnetic work, velocity/current derivatives, stress,
species work, metric initialization, induction solve and volume/face audits.
Methods and module aliases restore on success or failure. No numerical arrays,
equations, tolerances, acceptance checks or physical boundaries change.

The separate observation source archive includes all explicit dispatch modules
and the physical/guard closure. The physical source identity must be read with
this additional observation record. Fresh paths are mandatory. Existing records
and calculation sources remain immutable.

macOS `ru_maxrss` is cumulative process high-water RSS. A rise identifies an
interval establishing a new peak, not every live allocation or exclusive
function memory. Nested calls overlap. Resource/runtime metrics need not replay.
The existing 8 GiB RSS / 256 MiB free-disk guard applies; final observed values
must also satisfy both bounds. A resource interruption is not numerical failure.

This is a complete instantaneous physical workflow, not a trajectory, formation
experiment or stability/confinement claim. Its observations inform subsequent
workspace changes and guarded spatial refinement.
