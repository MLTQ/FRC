# Smooth shared finite-mass bulk transport

Reconstructs actual conservative rho, radial/axial momentum and angular momentum
with regular radial powers, plus separate pe, pi, current energy and magnetic
traces. Shared cubic normal momentum gives the mass flux. The same pointwise
mass flux multiplies velocity, angular momentum per mass and full specific bulk
enthalpy/kinetic/current energy. Face quadrature integrates nonlinear products;
mean-preserving transverse reconstruction retains the normal momentum average.
Axial face integration uses the annular r dr weight on every conserved component.

Maxwell/thermal stresses enter once. Staggered normal magnetic face averages
replace their reconstructed normal values; their transverse means are retained.
The radial hoop source uses the same pressure/rotation/magnetic reconstruction.
Full computed Poynting and drift stresses are added separately by the RHS.
The original energy partition and lower-order electromagnetic representation
are unchanged, so this is not a fourth-order claim for the entire model.

This is a central, smooth reference with no Rusanov or added physical diffusion.
No shock or nonlinear stability guarantee follows. Nonpositive density or species
traces reject; no floors, trace clipping, hidden fallback or physical atmosphere
are introduced. Current-energy interpolation can have signed tiny overshoots;
it is a numerical face trace, not a stored state or a claimed positive energy
quadrature. Production state current energy remains its positive physical form.

The axis/wall have zero normal mass, wall Btheta=0 and periodic axial faces.
Cache keys are immutable meshes and quadrature orders, and matrices are read-only.
Inputs remain unchanged. Compare contacts, all wave pressure rates, independent
manufactured forces, conservation and actual time evolution before adoption.
