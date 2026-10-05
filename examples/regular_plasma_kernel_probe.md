# Original kernel rejection capture

Runs the archived full coordinate-difference spectrum until its original strict
KernelInduction failure. A diagnostic wrapper captures the forcing, separated
kernel/remainder, already valid inner induction result, restored full E and both
curl evaluations, then re-raises the original failure. It never accepts that
rate or loosens a check, and always restores the original method on exit.
The latest capture also retains the exact fluid/magnetic input and component
spreads to identify the failing direction. The earlier capture is preserved in
regular_plasma_kernel_before_state; a pressure-only cause must not be inferred
from the generic exception without this input evidence.

Freezes the current coupled closure plus the old coordinate-spectrum source and
this diagnostic; archives array hashes and norm evidence without overwriting.
This reproducer supports diagnosis of frc-7mm.6. It is not a solver modification
or proof that a full nonlinear step should have been accepted.
