# Cylindrical original-platform evidence audit

Verifies frozen source and every raw-array hash. Rate and spectrum calculations
must reproduce arrays and metrics bitwise on the original platform. Trajectories
must reproduce initial state, final RHS, last accepted step, and every saved
physical history row from its actual state. A failed step must fail again
without mutating its accepted input. No failed case is promoted to complete.

Independent EOS/kinetic chain-rule energy rates are compared with production H
gradients. Initial/final regular physical integrations and 32/48-node,4x/8x
quadrature controls must replay exactly, with their difference explicitly saved.
The audit reports quadrature or accuracy limitations; bitwise replay alone is
not a physical accuracy certificate. Thresholds/convergence are assessed in
the verification report. Live checkpoints are excluded from terminal results.
