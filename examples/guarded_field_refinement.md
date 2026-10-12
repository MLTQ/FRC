# Single guarded full field-direction calculation

Dispatches the frozen direct-metric full signed workflow on F64/256, at the
requested field nr/nz, epsilon1e-4 and unchanged tolerance1e-10. Full actual
fluid/current/species terms, independent q16/20 audits, saved raw states/rates
and rejection behavior are inherited unchanged. No source or acceptance patch.

A separate control archive records own-process peak RSS (macOS bytes), minimum
free disk and whether a256MiB headroom guard interrupted execution. The guard
checks once per second; it is not a guarantee against transient exhaustion.
Its interruption cannot mark physical progress complete. The physical output
retains its original summary/schema/source identity. Only current-process
resources are read; no other process or global-memory interface is inspected.

Run original and repeat in separate processes. Optional --original invokes
the complete frozen all-array/non-runtime-metric verifier after the repeat.
A reproduced numerical rejection remains failed, even if replay verifies it.
Control runtime/resource measurements are not expected to replay identically.
Each destination/control path must be new. The control source closure is
separately frozen and checked after execution, including interruption.
