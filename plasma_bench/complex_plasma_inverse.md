# Complex shift for real two-stage collocation

Reuses the fixed physical block-preconditioner matrices and forms the Schur
inverse of I-h L for a complex h with positive real part. The nonlinear state
and equations stay real. For two-stage Gauss, h=dt(1/4+i sqrt(3)/12).

The real stage matrix has conjugate eigenvectors (i kappa,1) and
(-i kappa,1), kappa=2-sqrt(3). For real residuals r1,r2, one conjugate solve
with .5r2-.5i r1/kappa yields real stage corrections
(-2 kappa Im(w),2 Re(w)). No complex plasma state, extra damping or physical
term is introduced. Dense magnetic and fluid-to-magnetic coupling storage
scales with field DOFs; this initial preconditioner is for bounded verification,
not yet a scalable large formation-grid algorithm.
