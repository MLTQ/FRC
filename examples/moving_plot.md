# Scientific moving-plasma figures

`plot_errors` shows analytic sphere energy/force error, spurious translation
force, exact gas-expansion error, initialized-FRC drift and source-work time
refinement. Axes carry units and distinguish error from physical response.
Resolution/time axes label only the tested values; logarithmic minor labels
are suppressed to prevent collisions with those explicit ticks.

`field_panel` colors actual material Q1 polygons by P0 gas pressure and draws
flux contours from sampled r² times the actual Q1 a interpolant. The surrounding
gray region has no fluid. Contour triangles are display interpolation only,
not a solver/topology detector; collapsed triangles at the center are omitted
using their explicit signed meridional area.
`plot_fields` combines those fields with measured equilibrium drift and the
actual expanding surface's analytic radius error. No generated image or
prescribed animation substitutes for computed state.

Both PNGs must be visually inspected before reporting the study. Plot source
is hashed separately by the report and does not relabel measurement sources.
