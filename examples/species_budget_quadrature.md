# Common-grid quadrature comparison for species budgets

Reconstruct the actual q16-minus-q20 coefficient differences on an independent
common q20 physical grid. Compare every weak-budget term, each mass-diagnostic
rate and all refined-space rates. This avoids inferring field agreement from
nearly equal scalar errors. The full original entropy-rate norm supplies a
separate normalization for tiny or vanishing face terms. No norm floor.

`calculate(records)` expects the three complete verified input records, uses
their actual fluid/field dimensions and original physical species masses, and
returns a row-labelled integral table plus physical-difference fingerprints.
Columns are squared difference, squared term, squared full production rate.
The caller hashes and verifies inputs, archives source and performs a replay.
No solver, state change or tests. Tiny-term own-relative errors remain visible.
