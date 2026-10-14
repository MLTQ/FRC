# Physical samples of the deformed magnetic vacuum

Reconstructs the actual finite-space flux and magnetic field at physical points,
including the regular axis and the exact multipole exterior. Uses the same
moving Q1 geometry and complete polynomial flux/multipole trace maps as the
energy solve. No interpolation of a displayed quadrature cloud is used.

Inside the outer matching box, locate a physical exterior cell and invert its
bilinear coordinate map. Evaluate the reference flux and both derivatives there;
the inverse geometry Jacobian transforms derivatives into physical Br,Bz. On
the axis, the regular polynomial Bz limit includes radial stretch squared.
Beyond the fixed outer box, directly evaluate the saved decaying multipoles.
Points inside plasma or failed inversions reject rather than extrapolate.

This enables exterior-field visualization and the independent coil-linkage
comparison 2pi*sum(turns*psi_induced_at_coil) versus integral(J.Aunit). Those
finite-space quantities need not agree exactly; their discrepancy must not be
used to repair the source work or alter the physical currents. `geometry` may
be supplied to reuse the same already assembled vacuum state.
