# Open-boundary measurements

`static_measurements` compares asymmetric compact-current solutions to independent
loop quadrature. Spatial refinement uses interior nodal probes; modal refinement
uses analytic exterior magnetic fields, including an axis probe. Quadrature
refinement establishes the reference precision. `normal_derivative_jump` compares
second-order one-sided interior gradients to analytic exterior gradients at all
non-axis matching nodes, including corners. Its relative norm is unweighted.

`pulse_measurement` repeats the same complete two-coil pulse under selected
boundary, matching extent, mode count, mesh and CFL. Histories have the same 21
output times as the fixed-boundary report. `all_measurements` reruns both fixed
and multipole cases from current code; comparisons do not silently mix versions.
The full-space energy and beyond-grid energy are separately retained. Presentation
and provenance belong to `verify_open_boundary`.
