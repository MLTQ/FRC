# Repeated stationary-field measurements

The continuum reference has constant positive density/pressures, uniform Bz,
zero current and no thermal evolution. Eighteen controls span nr=8/16/32,
Bz=0.005/-0.005/0.5 T and dt=1e-8/1e-6 s for ten steps each. Field errors are
independently reconstructed from stored flux differences; the toroidal bound
uses the outer radius of each annulus. No weak-current threshold removes noise.

Retains every accepted residual, time, field maximum and conservative energy
change, plus initial/final physical states/current and hashes. Independent
physical quadrature checks total energy. Both species entropies must remain
bitwise unchanged, input arrays remain immutable, and archive values are reloaded
before acceptance. These controls test numerical equilibrium preservation under
the existing fixed-n, constrained-flow boundaries, not a free plasma interface.
