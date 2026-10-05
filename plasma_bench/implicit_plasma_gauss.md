# Full nonlinear fourth-order Gauss collocation

Two real stages satisfy Yi=U0+dt sum Aij F(Yj), with
A=[[1/4,1/4-sqrt(3)/6],[1/4+sqrt(3)/6,1/4]], b=[1/2,1/2]. All actual
fluid/current/species/magnetic terms use the same original complete RHS at each
stage. A complex shifted approximate inverse accelerates real chord iterations;
the approximate linear background never replaces an evolution term. A failed
chord line search tries a matrix-free Newton/GMRES proposal using the same full
residual and the same unchanged acceptance rules. Fallback counts and Krylov
histories are retained; exhausted Newton or line search rejects atomically.

Checks use actual rounded stage arrays in the original collocation residual,
the fixed physical active-coordinate infinity norm and the documented atol/rtol
rule. Returned endpoint rounding is checked separately. Actual midpoint/stage
and final physical states must remain finite and positive at points and common
quadrature. Invalid positive-state trials are damped; strict induction failures
propagate. No time retry, floor, filter, fitted field, heat or energy repair.

The linear stability function is (1+z/2+z²/12)/(1-z/2+z²/12), with fourth-order
time accuracy, unit modulus on imaginary eigenvalues and A stability. It is
not L stable: arbitrarily large steps need not damp physical fast decay or
resolve fast phase. Nonlinear energy is not exactly conserved. Full independent
wave/contact time/mesh and conservation evidence is required. Stage changes,
nonlinear/induction residuals and iteration histories are available for audit.
