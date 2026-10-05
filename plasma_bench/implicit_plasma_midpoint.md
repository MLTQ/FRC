# Full nonlinear implicit midpoint

Solves delta=dt F(U0+delta/2,B0+deltaB/2) with the actual complete coupled RHS.
The harmonic bias is fixed. A fixed uniform physical block preconditioner gives
chord iterations with a residual-decreasing line search; rejected positive-state
trial candidates are not accepted steps. There is no time-step retry, field
clamp, changed mass, heat, diffusion, projection or energy correction. If the
chord line search cannot decrease the residual, a matrix-free Newton/GMRES
direction is tried under the identical line-search and acceptance rules.
Fallback counts and linear histories are retained; linear failure is atomic.

Each residual evaluation retains all original strict induction checks. The
dimensionless infinity norm must satisfy atol+rtol max(|delta|inf,|dt F|inf),
using the documented fixed physical PlasmaCoordinates scales. Both midpoint
states and returned states must have positive point/sampled density and finite
physical EOS; returned rounded arrays are checked again in the original time
equation. Failure leaves caller state/field untouched. A failed inner induction
solve propagates; only invalid positive-state/EOS line-search trials are damped.

Implicit midpoint has second-order time accuracy and unit modulus for a purely
imaginary linear eigenvalue: (1+dt lambda/2)/(1-dt lambda/2). Large steps therefore
do not add linear damping, but do distort unresolved wave phase. Nonlinear
physical energy is not generally quadratic and is not exactly preserved by
this time method. Full wave/contact time refinement and conservation are required;
linear stability or nonlinear residual convergence alone proves neither accuracy
nor nonlinear long-time stability. Iteration, residual and induction diagnostics
are returned for each accepted step.
