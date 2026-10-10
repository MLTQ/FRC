# Streamed current and exact face audit

Uses the original independent current, strong magnetic curl and exact one-sided
trace formulas from radial_mode_interfaces. Volume current/gradient/target
arrays are reconstructed in radial batches; all quadrature points remain.
NormAccumulator and per-row math.fsum change only diagnostic reduction order.
The strong-curl comparison still excludes distributional interface terms.

Exact radial/axial left/right traces, the periodic seam, original analytic
denominators, jump maxima, per-face squared integrals and fingerprints are
unchanged. These face arrays are smaller and remain evaluated together.
No averaging, shifted face sampling, physical correction or numerical filter.

Verify through complete signed original-versus-streamed physical workflows,
quadrature controls and fresh full replay. No unit-test suite is needed.
This audits initial spatial representation; it does not establish complete
local energy evolution, vacuum coupling or particle confinement.
