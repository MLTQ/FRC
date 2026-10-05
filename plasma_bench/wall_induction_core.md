# Positive induction on admissible trial vectors

Parameterizes every magnetic-rate Krylov vector with fixed zero axis/wall psi
and a radial null-space basis of the exact toroidal wall trace. This avoids the
small inadmissible roundoff leakage observed when solving in the larger space
and hoping the exact solution's constraint survives iterative arithmetic.
It never projects a computed accepted rate or state after the solve.

With lift T normalized by its magnetic energy diagonal, solves
T.T [H+Aw.T I(n) Aw] T y = T.T H Fw M^-1 f, Bt=T y,
Aw=-M^-1 Fw.T H. The cached preconditioner is the exact restricted magnetic
Hessian; current inertia always uses the current physical density/masses.
The complete weak Ohm/Ampere and *full*, unreduced Fw Faraday physical residual
must still pass the unchanged tolerance. A small reduced residual alone does
not accept a solve. Inner correction/refinement budgets match the original
positive induction core, with atomic failure and no state mutation.

StrongWallInduction supplies separate analytic curl-kernel extraction/restoration
and original raw-incidence defect reporting. This core is for bounded wall
verification only; it is not vacuum or physical coil/end-boundary modeling.
