# Source, physical and trajectory audit

Checks frozen/live source closure, archived source bytes and raw data hashes.
Replays every rate array exactly and independently recomputes local physical
thermal derivatives at quadrature 18. A 1e-7 agreement ceiling is an audit guard,
not an accepted continuum error or a replacement for measured values.

Trajectory audit replays initialization and the last accepted RK4 step, checks
complete positive/finite accepted histories and final work increment, and
independently integrates all five saved states at order 18. It does not claim
to replay every preceding step. Continuum characteristic norms and independent
time-refinement tests supply accuracy evidence beyond algebraic conservation.
