# Refined induction checks

An independently initialized localized toroidal magnetic perturbation reproduces
the old original-tolerance failure. The opt-in refined assembly must pass every
unchanged original residual with a tighter inner solve and untouched input.
A successful smooth wave must return all numerical arrays and checked residual
bitwise unchanged. A deliberately insufficient iteration budget must exhaust
the refinement attempts and reject without changing the accepted state.
Full coordinate-Jacobian and nonlinear trajectory evidence is a separate study.
