# Physical thermal source, equation and trajectory audit

Checks live/frozen transitive source, archive and raw-array identities. Rate
inputs and both species' scalar solve/pressure fields are reconstructed exactly;
original full Ohm/Ampere/Faraday equations are checked after restoring potentials.
Full rate-relative Faraday error remains separate from backward error, including
absolute error when the target rate is zero. Independent quadrature14 replays
local physical density/entropy/current/magnetic energy derivatives and all norms.

Transport audits reproduce initialization and the **last accepted RK4 step**,
check time-history completeness and positive finite accepted pressures, and
independently recompute every saved snapshot's physical norms/integrals. They
do not claim a bitwise replay of every earlier time step; independent continuum
targets and matrix-exponential tests supply time/space accuracy evidence.

The1e-7 local-integral agreement ceiling is a quadrature/representation audit,
not an allowed continuum error or relaxed original solver tolerance. Every
measured defect and hash is retained rather than only a pass flag.
