# Uniform physical block preconditioner

Builds a fixed-background approximate coupled Jacobian to accelerate a full
nonlinear residual solve. Sparse regular fluid acoustics/advection are coupled
to a dense compatible Hall/inertial magnetic block and the exact reciprocal
Lorentz/ideal work adjoint. The magnetic metric includes physical current
inertia. Background density and pressures come from the integrator's documented
initial means; this background is never substituted into evolution.

Uniform axial magnetic translation uses a centered coefficient derivative as
an explicit approximation in the preconditioner. The actual nonlinear residual
still uses complete CurrentTransport and both species. Matrix construction calls
the strict refined induction solver. No artificial mass, energy correction,
viscosity or dropped term is introduced into the evolution equations.

`factor(h)` forms the sparse fluid inverse and dense magnetic Schur complement
for I-h A. `solve` applies that coupled block inverse. Cached factors depend on
h only; a new preconditioner is required for another initial geometry/background.
The approximation must be verified against separate linearization at rest and
nonlinear implicit solves; preconditioner convergence itself is not time accuracy.
