# Full captured-state curl evaluation audit

Reconstructs the original compatible geometry, wall projection and exact
DirectMetricWallInduction lift expressions without constructing an unused
preconditioner or performing another plasma solve. Each of the four archived
iterates must reproduce RHS, magnetic rate, dual residual and electric-curl
distribution vector bitwise before alternative evaluations are interpreted.

Measures assembled, sparse/dense factored and subtraction-first actions on
every complete e0, ei, e and independently verified addition-low array. At the
final iterate, independent exact Decimal raw products/sums provide a reference.
Every action and residual difference is saved. Physical Faraday uses the
unchanged original magnetic Hessian and actual magnetic-rate norm.

Alternative dual residuals use the original RHS normalization for comparison;
they are not accepted solver residuals. This workflow changes no solver or
state. Sparse projected-matrix assembly and each evaluation order remain
distinguished. Exact raw products do not certify an exact projected result.

The source closure, input hashes, immutable archive and complete separate
replay are required. A pre-existing output directory rejects. The unchanged
filled-periodic wall fixture establishes neither vacuum formation nor physical
component convergence. No unit tests, diffusion, filters or energy repairs.
