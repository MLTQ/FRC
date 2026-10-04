# Original current-advection and drift-stress verification

Initializes the unchanged physical regular-current fixture and evaluates each
new transport term separately against independent covariant derivatives of
analytic fields. Actual Ampere current and actual fluid continuity are used in
the candidate. Field/current projection errors, fluid reconstruction errors,
pointwise continuity defects and normal-current jumps are explicitly measured.

Face-omission controls preserve skew volume energy identities but expose local
force/electric accuracy failures. The complete ideal/Hall/inertial-transport
instantaneous solve checks whether physical density exchange is now paired;
species pressure/current entropy transport and time evolution remain absent.
No artificial heat or force correction is allowed to close the total.

## Historical pause boundary (2026-09-30)

Candidate only. `run_case()` has focused-test coverage; `main()` is unfinished
because `current_transport_records` has not been written. No frozen 16-case
study, independent audit, or full regression has run for this candidate. Resume
with beads frc-7mm.3 and docs/NEXT_SESSION.md before invoking the study CLI.

## Resumed verification (2026-10-04)

Records/provenance and independent audit are implemented in companion modules.
The CLI completed all sixteen cases under source 2f682eec; arrays and metrics
replay bitwise. Full regression has 408 passing tests. Independent continuous
force/energy diagnostics and measured limits are in CURRENT_TRANSPORT_VERIFICATION.
The paragraph above records the earlier pause, not the current CLI status.
