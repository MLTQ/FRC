# Independent mixed-flow audit

Source and raw hashes, initial state, final original RHS, last accepted RK state
and final source-work/source-total increment replay bitwise. Failed steps must
fail again with unchanged accepted input. Every saved physical history is
recomputed from its actual state; the accumulated source-work ledger before
the final step is retained evidence, not claimed independently fully replayed.

Every accepted snapshot is reconstructed on a common physical grid and a
higher32-element/40-Gauss-node/512-axial grid. Continuous energy/angular history,
grid-versus-continuous differences, sampled positivity and quadrature differences
are retained. Exact manufactured field errors replay independently. The physical
kinetic/EOS energy rate is checked against production H-gradient work.
An unforced trajectory has no manufactured final target. Raw convergence and
time-control comparisons are separate from the audit's replay assertions.
