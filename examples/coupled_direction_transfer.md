# Independent differential and finite-state transfer audit

Verify nested conservative subdivision of both full archived tangent directions
onto F16/64 and F32/64 at fixed B32/64. Original normalization, magnetic/current
coefficients, physical masses, axial nodes and stationary base remain unchanged.
The radially constant base is checked before exact replication.

At q16/20 independently reconstruct density, all velocities, both entropies,
both pressures, B and current. Differential transfer must agree within 1e-9
relative for every nonzero physical component; zero targets retain absolute
errors. This is a transfer diagnostic bound, not a solver residual tolerance.
Every owning-parent conservative face is retained bitwise, including the
inactive-axis and normal-wall contracts. B/current bytes and samples match.

For each sign at h=1e-4,1e-5,1e-6, independently compare the complete finite
physical states too. The nonlinear nodal q/n velocity interpolation can change
under conservative subdivision. Those differences are normalized by h times
the original differential component, retain absolute errors, and are not hidden
as exact primitive transfer or replaced by a projection.

Save both fine contexts, every metric/fingerprint/input hash and complete source
closure. No coupled RHS is evaluated until the complete transfer protocol passes
and a fresh replay reproduces all raw arrays and metrics. No unit tests, evolved
state changes, numerical diffusion, density floors or mass alterations.
