# Checkpointed common-radial-field trajectories

Fourteen unforced verification fixtures retain the contact/full-transit,
pre-shock nonlinear acoustic wave and rotating equilibrium physical inputs.
Contact steps are half the historical unsuccessful coarse step; nr8 quarter-step
and reverse-flow controls remain explicit. Acoustic nr16 half-step, reverse
and higher quadrature controls are included. These are not formation experiments.

Actual face radial momentum, cell density/axial/angular momenta and both physical
entropies advance without forcing or repairs. Separate histories retain physical
energy, conservative energy account, RK4-integrated spatial energy defect,
positivity, mass/momenta, both entropies and solver residuals. Five snapshots,
initial rates and the last accepted step are retained. Failed steps preserve
accepted arrays and the error; failure never implies completed duration.

At roughly 5% intervals, accepted state, complete history, snapshots and energy
ledger are written with immutable array filenames and an atomic metadata pointer.
`--resume` validates source/raw hashes and bitwise initial state/rate, then
continues only incomplete cases from their actual accepted step. It does not
retry a terminal failed case. Resume requires the original numerical platform.
The latest accepted step between checkpoints can still be lost on interruption.
