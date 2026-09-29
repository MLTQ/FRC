# Physical magnetic-rate verification after inertia

Solves the actual static discrete Ohm/Ampere/Faraday system with physical inertia,
separately for full original forcing and only its nonlinear convective term.
The auxiliary tau rescales linear unknowns; it is not time integration or an
effective mass. Both RHS share one checked linear factorization. The result is
physical gdot, unlike the earlier intermediate B* forcing measurement.

Independent Fourier/Chebyshev continuum solves at radial orders 64/128/256 and
256 axial points supply reference means. Actual collocation residual, analytic
axis-probe sensitivity and inter-resolution reference differences are recorded.
All five production meshes retain volume/axis/maximum errors and lossless fields.
Engine and every extra reference/driver source are frozen and hashed. This is
an instantaneous consistency study, not nonlinear moving-plasma evolution.
Each continuum and discrete archive is hashed before inclusion in the final
collection; a syntax or solver failure cannot generate a completed report.

`--method schur` selects the algebraically equivalent positive magnetic-space
elimination in current_rate_schur, avoiding expensive mixed LU at fine meshes.
The method and its additional source hash are recorded explicitly. The default
remains `mixed`; use a separate output directory for the alternate study.
