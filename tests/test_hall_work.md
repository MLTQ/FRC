# Experimental Hall spatial work checks

An independent finite difference of the declared magnetic energy validates its
gradient and work conjugate current. It uses a resolvable centered perturbation
for the exactly quadratic energy;
there is no finite-difference truncation error to trade for cancellation.
Nonlinear variable-density snapshots verify global Hall cancellation, both
Faraday components, CT, conducting and
periodic boundary fluxes, and strictly positive carrier preconditions. Local
Poynting mismatch remains explicit; it is never added as physical heat.

Physical-volume error against independently differentiated analytic fields
refines globally. The test deliberately makes no axis accuracy or nonlinear
evolution claim: the study records a persistent first-annulus gdot defect.
The experimental operator is not connected to ToroidalMHD. Fixing that defect,
local transfer, stable time stepping and full wave evolution remain required.
