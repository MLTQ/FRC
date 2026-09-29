# Nonlinear toroidal measurement records

Runs all declared cases and preserves independent continuum comparisons, complete
conservative state/flux snapshots, initial RHS, physical coordinates/volumes and
accepted-step diagnostics. Histories are lossless numeric columns to avoid huge
repeated object overhead. Records include exact source/config/archive hashes and
dependency versions. Source identity must stay fixed throughout every run.

Errors compare declared cell averages with independently constructed references.
Rotation and Btheta use physical r² volume weighting; pressure/density use volume.
Absolute errors remain available when a reference field vanishes. Decay records
include continuum Bessel heat, near-axis g error and global toroidal field error;
its nonlinear feedback is asymptotically small, not identically absent. Exact
nonlinear contacts use true continuum E, retaining subcell recovery discrepancy.
No incomplete case is marked complete; reporting verifies every identity.

Physical tensor-Gauss norms of the actual limited reconstruction against true
continuum point fields are also recorded, with separate orders 6/10. They
include subcell error and are explicitly separate from cell-average norms.
