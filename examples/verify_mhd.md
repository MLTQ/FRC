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
first-axis-cell absolute acceleration error, keeping local axis accuracy
visible alongside global conservation and volume norms.
Refinement plots label actual cell-count ticks and suppress crowded minor labels.
The default measurements now use second-order reconstruction. The companion
mhd_accuracy module isolates smooth and limiter-switching temporal cases and
reruns first-order comparisons at identical mesh/parameters. Reports explicitly
distinguish the special straight-field equilibrium from the independently
verified curved equilibrium and its measured drift. Low-pressure results retain
disabled-reference failures, limiter activity, stronger stress cases and cold
rarefaction errors. Positive traces are distinguished from full-update accuracy.
