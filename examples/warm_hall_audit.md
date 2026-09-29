# Complete warm Hall audit

Requires all 58 exact configurations, current transitive source identity,
matching per-case/assembled JSON, archive SHA256, finite arrays and coordinates.
Quarter-cell sums independently check all stored dual volumes with a declared
radial-subtraction roundoff allowance proportional to nr machine epsilon.
Ion-pressure replay uses a 32-epsilon allowance scaled by the separate energy
terms being subtracted, so nearly depleted positive ion energy is not tested
against an ill-conditioned relative-to-small-residual roundoff criterion.
Replays every compact history endpoint/count, positive species minima, monotonic
time, actual step caps, nonlinear tolerances, species sum and conservation drift.
Physical logarithmic entropy is measured separately from the linear proxy.

Independent field audit covers full moving endpoints. Isolated pressure runs
replay electron, ion and toroidal magnetic energy and entropy from saved fields.
Standalone diffusion compares physically integrated field loss with recorded
positive heat. The independent analytic face-flux divergence is compared with
analytic local energy work, including barotropic cancellation, without calling
the production energy operator. Exact finite-volume flux replay is a separate
archive integrity check. Failing assertions never repair or relabel records.

An explicit configuration list can audit a supplemental successful subset
directly from its individual records. The supplemental report must separately
verify complete protocol coverage including rejected cases; this mode does not
create or publish an assembled complete collection.
