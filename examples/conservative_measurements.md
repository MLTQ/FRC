# Conservative vacuum-coupling measurement suite

Records original standing shocks, first-order/limiter/interface/time controls,
axial-boundary sensitivity, identical prescribed-coil pulses and matching-surface
controls, plus independent smooth advection/diffusion refinement. The original
sharp shock's remote boundary data are retained; no pressure buffer or easier
replacement initial state is substituted. Disabled-method failures remain
records with accepted time and diagnostics, not omitted data points.

Shock errors use the previous study's exact cylindrical overlap-volume window
and independently averaged reference. Source histories have 21 output times,
including absolute response residual and error normalized to actual input work.
Fluid energy-flux conservation, independent response error, interface numerical
heat, physical heat and flux-limiter activity are separate diagnostics.

Each case is saved when complete; explicit reuse requires identical package,
measurement source and environment. Complete parameters, source waveforms,
initial/final diagnostics, history and radial shock profiles are persisted.
Failures must not be treated as final-time convergence points in the report.
Sparse per-cell accepted limiter maps record affected and troubled Euler-stage
counts, grid centers, maximum distance from a gas boundary and any direct
limiting in the independent shock's reference window. Empty entries mean zero,
not missing samples. Aggregate counts remain in every history diagnostic.

An optional `progress` observer receives accepted-state summaries at setup,
each requested output and at most every 30 seconds during advancement. It does
not change steps or physical parameters. These records are explicitly marked
`status=running, completed=False`; they are progress evidence, not restart files
or final convergence points. The compatible-energy driver persists them in
separate progress files until the final complete record is saved.
