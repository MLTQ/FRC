# Resource guard for the complete weak entropy audit

Explicit opt-in counterpart of `guarded_spatial_parts`, dispatching only
`entropy` to `species_interface_audit.run`. Both calculation and monitor are
included in the frozen source closure before execution. The monitor checks
own-process macOS cumulative RSS and free disk each second, with final extrema.
Unchanged limits are 8 GiB RSS and 256 MiB free disk. A resource interruption
preserves partial evidence and cannot be called a completed scientific audit.

`--original` performs the complete non-runtime-metric and raw-NPZ replay verifier
after a fresh process rerun. All raw hashes must agree, not only array values.
The numerical identity reporting threshold is distinct from physical solver
tolerances; this workflow runs neither induction nor time evolution.

Run one heavy worker at a time. Collect actual process exits before calling the
work complete. The monitor's status file alone does not establish liveness.
No unit tests, source changes, retry, physical state projection or repair.
