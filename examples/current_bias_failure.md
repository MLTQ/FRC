# Uniform-bias precision observation

Uniform n, pe, pi and axial B is an exact stationary continuum state with J=0.
Two meshes and two requested time steps preserve actual solver success/rejection,
tiny represented weak-current values, accepted time, unchanged-input verification,
physical inputs and lossless endpoint arrays. Rejection is not silently treated
as successful zero evolution. This exposes a precision/residual acceptance gate
that must be addressed before using the component for general initial states.

The repeated current difference in the returned-state inertia residual is a
candidate cancellation mechanism, not a proven attribution. A fix must retain
honest finite state/residual accounting and strong-field small-perturbation
controls without arbitrary heat, mass adjustment, pressure floor or state reset.
