# Finite-cylinder inertial induction

Uses the same strict VacuumInduction/SplitInverseWallCore equation and bounds.
The lift includes every nonaxis flux node, including both ends and the side;
toroidal coefficients use both radial and axial trace-null bases. The magnetic
metric includes eliminated exterior energy. Current inertia uses actual density.
Uniform density is only the Cholesky preconditioner. Kernel splitting preserves
axially constant Er and radially constant Ez on the finite nodal/edge dimensions.
No periodic trace, state repair, numerical diffusion or tolerance change.

Seven full24us coil pulses now exercise this coupled path. Independent limits,
spatial/time comparisons, transport flags and remaining exterior-resolution
limits are recorded in docs/FINITE_CYLINDER_RESULTS.md.
