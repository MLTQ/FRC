# Independent archival and physical replay

Requires full ordered protocol, unchanged source, exact individual/collection
JSON, NPZ hashes, finite arrays, mesh coordinates and complete accepted times.
Success and rejection have distinct time/state contracts.
Accepted final time permits only two floating-point eps relative rounding from
the declared integer number of equal steps, with no early-finish tolerance.
Endpoint physical energies are integrated independently at orders 8/12; separate conservative
cell energy and inferred ion pressure are replayed, including rejected trials.
Density/current kinetic energy use positive Q1 density, not a cell shortcut.
Pressure-defect replay records an absolute roundoff allowance from positive
energy accumulation and subtraction: gamma_k times the sum of conservative and
component energy magnitudes, in pressure units, with
`k=2*(8^2+12^2+production_order^2)+64`. This conservative operation budget covers
the compared quadrature accumulations and scalar arithmetic. It extends the
original 1e-10 Pa audit allowance when magnetic energy is much larger than the
pressure defect; it is uncertainty in replay, not an accuracy or positivity bound.

Independent dual-volume supports replay both species entropy and pressure.
Initial continuum samples, full-electric Faraday differences, boundary powers,
exact-wave norms and face-integral quadrature controls are checked separately.
The full-electric rate replay also records absolute discrepancy and the storage
roundoff bound `8 eps ||(|final|+|mid|)/dt||`, covering the state addition,
midpoint, subtraction and division. Relative discrepancy must pass the original
2e-7 threshold or fit this explicit bound. Weak waves on large backgrounds can
exceed that relative threshold even when their state matches Faraday to storage
precision; this is recorded, not claimed as improved small-signal resolution.
The new stored-state residual history, when present, must pass its tolerance.
All-history global drifts are retained alongside local energy defects. The
smooth-protocol quadrature tolerance is not a claimed arbitrary-density bound.
This audit is read-only and neither masks rejections nor thresholds physical
wave errors into a generalized accuracy claim.
