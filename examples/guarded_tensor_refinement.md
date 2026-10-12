# Guarded fine-grid tensor candidate

Runs the complete opt-in tensor-mass signed workflow on F64/256 with requested
field nr/nz, epsilon1e-4 and unchanged tolerance1e-10. Original and replay use
separate processes; --original verifies all arrays and non-runtime metrics.
The separate control source archive records resource and runtime values.

Preflight and a once-per-second monitor enforce256MiB free disk and8GiB
own-process peak-RSS budgets. These are execution budgets, not altered plasma
physics. If breached, request SIGINT, then exit only this worker with code130
after5s if native code does not yield. The control records the interruption
before forced exit; physical summaries remain explicitly incomplete. No other
process is inspected or terminated. This fixes the old guard's inability to
finish cancellation while native factorization ignored a pending SIGINT.

Resource checking cannot guarantee protection against changes between samples.
Runtime/resource metrics are not expected to replay bitwise. A resource event,
strict numerical rejection and successful physical workflow remain distinct.
All old sources and interrupted evidence stay unchanged. No unit tests.
