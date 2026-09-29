# Compatible Hall scientific figures

Three standalone Matplotlib figures show the corrected and historical regular
axis errors, nonlinear spatial/current/local-energy errors, actual three-branch
wave propagation, fixed-grid time convergence, restricted nonlinear mesh
differences and independent conducting Bessel field/heat errors. The caller
must audit records before plotting. All axes state their norm and resolution.

Initial physical magnetic representation errors are dotted on the wave panel,
so the first-order field representation is not mistaken for accumulated time
error. Slow and other branches have different stated durations. Coefficient
errors, physical field errors and restricted conservative-state comparisons
are separately labeled. Figures do not imply formation or thermal closure.

`figures` returns output paths for the report's separate SHA256 manifest.
