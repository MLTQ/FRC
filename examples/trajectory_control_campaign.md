# Serial complete-trajectory controls

`run(destination, plan_path)` runs frozen guarded_trajectory_controls jobs in
fresh child processes, one at a time. Each plan entry supplies a simple name,
physical/control destinations, four mesh counts, signed epsilon, step count and
optional original trajectory for full replay. Duration and original numerical
bounds remain fixed by the guarded driver. Plan and source hashes must remain
unchanged throughout. Existing destinations are refused by the child.

Before each child, require256MiB free disk plus a declared planning reserve.
For the existing degree-five, two-stage, equal-axial-count protocol, derive raw
array counts from the original full state, packed magnetic field and reduced
stage coordinates. Reserve1.1 times the entire uncompressed trajectory payload
plus8MiB for metadata/archives. This is a conservative planning estimate, not
a proven upper bound; each child still enforces the original continuous disk/
RSS guards and final extrema. Other axial-grid relationships reject until a
separate storage derivation is provided. Git storage is not in this estimate.

A resource preflight stop records an unstarted job and aborts the remaining
queue. It is not a physical rejection. Nonzero child exits also stop the queue
and preserve exact partial evidence. A terminal numerical trajectory failure
can advance to its separately declared full replay; it is never marked as a
successful100us trajectory. all_jobs_terminal and all_full_duration_complete
are distinct. Only the former can hold after reproduced numerical failures.

The parent records actual child PIDs, return codes, log hashes and complete
terminal guard summaries. Interruption signals and collects only its own child.
Parent status summaries are written by atomic replacement.
There is no automatic time retry or relaxed bound. After complete verification,
the existing hash-checked sharing helper shares identical completed replay
NPZs while preserving all paths/bytes and marking them immutable.

The parent is orchestration, not an independent physics audit. Opposite-sign
physical/conservation audits, same-grid arithmetic controls, mesh/time analyses
and their full replays remain necessary. Do not run alongside another large
solver/audit. No prior source, physics, state, numerical diffusion, floor,
atmosphere, mass or energy repair is changed. No unit tests.
