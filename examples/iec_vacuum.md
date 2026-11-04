# iec_vacuum.py

## Purpose
Run the first open-cage 3D vacuum calculation for frc-cac. It varies numerical
resolution for a specified geometry, not plasma operating parameters. No
physical plasma, electron source, pressure, drive frequency or fusion model
is present; the historical air discharge is not imposed as a design target.

## Components
- `spheres`: an independent concentric-conductor Gauss-law potential, field
  and capacitance comparison at increasing harmonic-basis resolutions.
  This is a verification fixture, not a cage or formation experiment.
- `finite_line`: independent free-space finite-line potential and bisector/
  axial fields, exercising the analytic line kernel and the central image limit.
- `interior_field`: separate radial force samples and decreasing-difference
  divergence residual in empty space. Centre E=0 by itself says nothing about
  trap curvature; radial-force signs are retained without angular averaging.
- `main`: constructs the whole reference wire cage and measures fitted and
  separate off-fit surface errors, chamber-potential error, image quadrature
  convergence, centre field/potential and capacitance. Raw coefficients,
  source hashes and local copies are retained at entry with the results.
  All resolutions use the same off-fit surface point set, and coefficients
  from every refinement are saved. Changing assessment samples must not be
  mistaken for convergence of the boundary approximation.
- `plot`: depicts real wire edges, a vacuum potential plane and separate
  radial rays. It does not average the field or assume a plasma potential.

Usage:
`OPENBLAS_NUM_THREADS=1 .venv/bin/python -m examples.iec_vacuum OUT`
with optional `--panels 2 4 6 --diameter .15 --wire .002 --wall .2`.
The output directory must be new; earlier measurements are never overwritten.
`--wall` is the radius of the assumed grounded sphere, not its diameter.
Metal exclusion uses the rounded-wire union. Chamber/feedthrough/CAD details
are not known, so results are for an explicitly labelled reference geometry.

## Interpretation
The field satisfies vacuum Laplace outside interior mathematical sources.
Finite source and surface representation limits the equipotential accuracy;
use off-fit errors and actual refinement rather than the fit alone. A nearly
flat centre in the regular empty cage is not evidence of plasma trapping.
Kinetic electrons and ions, actual electrode supply, collisions and loss/work
budgets are still required before scanning for useful POPS parameters.

Measured results: `outputs/iec_vacuum_reference_v4/report.json`. The independent
sphere comparison's potential max error decreases 6.89 -> 2.16 -> 0.519 ->
0.0545% over 12/24/48/96 sources; field L2 13.5 -> 4.46 -> 1.02 -> 0.142%.
Capacitance relative error is 9.7e-7 at 96 sources. Finite-line analytic fields
agree to 2.3e-16. The empty-space field divergence diagnostic decreases
2.13e-5 -> 5.31e-6 when its difference spacing halves.
The cage's local surface-potential approximation is still limited to roughly
1.4% worst case; its error is not inferred from the sphere comparison.
Earlier coarse fits had small fitting errors but large off-fit errors and
are retained as preliminary evidence, not used to claim electrode accuracy.
