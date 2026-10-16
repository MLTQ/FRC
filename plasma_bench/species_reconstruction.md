# Conservative species face reconstruction

Linear reconstruction acts on the conservative vector (number, three scaled
momenta, total translational energy), so its two endpoint values average exactly
to the original cell state. Independently reconstructing primitive pressure and
velocity did not retain this kinetic-energy consistency and failed in strong
moving rarefactions even with physical endpoint pressure.

The nonuniform-grid componentwise minmod gradient uses one-sided endpoint
gradients. A single factor theta in [0,1] multiplies all five conservative slope
components for each cell/species. If both traces are already admissible theta=1.
Otherwise bisection on the convex Euler admissible set finds a nonnegative
number/thermal-energy pair. Exactly empty species remain exactly empty; no
positive density floor is introduced. Only face slopes change, never cell
numbers, momenta or energies. Limiting adds disclosed numerical diffusion.

This uses the conservative slope-scaling principle of
[Zhang and Shu (2010)](https://www.sciencedirect.com/science/article/pii/S0021999110004535).
It does not claim that their complete DG theorem applies to this moving-grid
RK4/electric-source discretization. Actual intermediate and final states must
still be valid; rejected states are preserved without repair. The returned
theta permits callers to report how often reconstruction was restricted.
