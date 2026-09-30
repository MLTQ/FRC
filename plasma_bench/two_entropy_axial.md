# Periodic axial verification of the joint physical flux

Uniform nodes store point values, never actual cell means. The skew periodic
centered derivative of order2/4/6 is applied by symmetric two-point flux
differencing: U_t,i=−2 sum_j D_ij F(U_i,U_j). This preserves all six state
integrals and physical H under the nodal dx integration norm. The corresponding
local energy-flux divergence is retained separately from w dot U_t.

This high-order composition follows the general centered flux-combination
framework described by [Klein and Öffner](https://arxiv.org/abs/2211.01171).
The joint two-entropy flux is separately derived in two_entropy_flux.md.

RK4 advances every variable without a temperature or energy reset. Its finite-
time energy defect is not hidden by the spatial identity. Positivity is checked
at every stage and final state; failure leaves the input unchanged. No limiter,
viscosity, entropy production or positivity theorem is added. Periodic axial
tests are a step toward joint transport, not a substitute for cylindrical,
magnetic/current, free-edge/exact-vacuum and formation requirements.

The energy theorem concerns the declared nodal quadrature. Independent denser
integration of a reconstructed physical state must assess aliasing and cannot
be replaced by this exact discrete invariant. Constant-pressure contacts and
full nonlinear waves also require independent continuum accuracy checks.
