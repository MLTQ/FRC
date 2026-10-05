# Modal preconditioner and actual evolution controls

All Fourier modes of random physical fields, including zero/Nyquist, are
compared with the separately constructed full-grid uniform physical blocks at
rest and with axial motion. Complex modal corrections must satisfy the full
real coupled two-stage Gauss equations, not just the block-factor implementation.

Centered differences of the actual complete nonlinear uniform-rest RHS must
converge quadratically to the separately assembled modal parent operator. The
1e-5/1e-6 pair resolves the central-difference limit; a preliminary 1e-4/1e-5/1e-6
scan measured 2.16911e-6/2.16913e-8/2.17059e-10, confirming second order rather
than requiring a coarse finite-difference sample to meet a tighter accuracy.
Both ordinary and secant-accelerated actual
wave stages/endpoints are independently reevaluated against full nonlinear
equations and boundary constraints. Deliberate contact iteration exhaustion
must preserve both input fluid and magnetic arrays exactly.

These tests verify equivalence and solves, not long nonuniform trajectory
accuracy or formation. Independent growing-mode/time/mesh controls remain needed.
