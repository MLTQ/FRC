# Coupled Hall-fluid measurement

Evolves the complete conservative HallMHD state, saving every accepted diagnostic
history and full initial/final arrays. Sources/configuration must remain fixed.
Wave initial states use true continuum averages; reconstructed velocity,
density/pressure and compatible magnetic fields are compared at actual physical
quadrature points. Magnetic and thermodynamic errors are normalized by the
wave perturbation, not a large constant background. Initial spatial errors and
final errors are separate. Conservative g cell errors remain a distinct metric.

Wave references are linear, so finite-amplitude residuals are explicit. The
nonlinear manufactured initial conditions have no claimed exact trajectory;
their separate refinement and time controls compare actual conservative states.
Histories preserve positive pressures/densities, source-free mass/energy/torque,
CT, physical heat and nonlinear/linear solver residuals. No failed or partial
evolution is saved as completed. A zero-mean angular-momentum total is not used
as a misleading relative-error denominator.
