# Matrix-free nonlinear correction proposal

When the fixed-background chord direction cannot decrease the actual full
implicit residual, a centered finite difference of that same residual supplies
Jacobian products to preconditioned GMRES. The physical coordinate perturbation
is cbrt(machine epsilon) times (1+max absolute stage increment). The existing
uniform physical block inverse is only a preconditioner. The linear proposal
uses relative tolerance 1e-4, up to five restart cycles of 30 directions.

This does not change the nonlinear or induction acceptance bounds. Every
Jacobian evaluation uses the complete original strict plasma equations; the
candidate must still pass the original positive-state, residual-decreasing
line search and final collocation equation. Linear failure rejects atomically.
No density, pressure, current or field component is discarded or repaired.
Finite-difference truncation and roundoff affect the proposed direction; they
cannot authorize an otherwise rejected time step. Record Krylov counts and
residual history so fallback cost and failure remain visible.
