# Full Radau plasma trajectory

Uses ImplicitPlasmaRadau with the independent exact contact or finite-mass
wave initialization, through a complete requested transit/period. Records all
physical field errors, grid conservation, residuals/limits, iteration and Newton
counts and accepted progress arrays using the same conventions as the original
implicit trajectory driver. Every actual stage uses the full physical RHS.
No rejected trial is recorded as a completed time step.

This separate driver preserves the frozen Gauss/midpoint source closure and
its failed-contact evidence. Radau supplies declared temporal damping, which
must converge under time refinement; it is not physical transport. Independent
continuous quadrature and full replay live in radau_plasma_run's audit action.
