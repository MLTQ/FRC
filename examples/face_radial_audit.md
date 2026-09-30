# Independent saved radial-face audit

Requires the exact source identity and raw hashes. Recreates initial state/rate,
every rate-study array and every terminal trajectory's last accepted RK4 step
bitwise on the original numerical platform. Failed cases also replay their
actual rejected step with unchanged accepted arrays.

A separate order18 root-product/Legendre sampler integrates every snapshot and
final physical state. Higher-quadrature local energy agreement and continuum
force/thermal norms are independent of the production matrices. Failed-state
quadrature warnings are reported explicitly. History totals replay separately.
Physical energy drift is separated into conservative-account drift, accumulated
spatial defect and remaining finite-step chain-rule error; none corrects the
evolved temperatures. Running checkpoints are not terminal verification records.
