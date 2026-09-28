# Quantified dynamic verification
`verify` runs independent fast-wave, 2-D advection, resistive-decay/heating and
Sod benchmarks over several grids, retaining parameters, source fingerprint,
dependency versions, measured errors, conservation diagnostics and plot arrays.
`plot` renders a scientific comparison with analytic solutions. `report` writes
the measured findings without claiming physical formation validation. Finite
resolution and numerical heating discrepancies remain visible. Paths are
explicit under outputs/mhd_verification and docs/MHD_VERIFICATION.md.
A separate fixed-grid time refinement measures RK2 accuracy against a finer
temporal reference, without confusing spatial and temporal truncation error.
A known static magnetic-pressure equilibrium measures spurious flow caused by
numerical diffusion across refinements; the report does not call it an FRC.
Analytic magnetic-force comparison includes both volume-weighted error and the
first-axis-cell absolute acceleration error, exposing weaker pointwise axis
accuracy rather than hiding it in a global conservation number.
Refinement plots label actual cell-count ticks and suppress crowded minor labels.
