# Constrained-current energy flux

Reconstructs full E from the reduced induction field plus the pressure potential
and -delta mu grad(|w|²)/(2e). The last term is the kinetic gradient required by
the vector identity for convective inertia; zero curl does not remove its
Poynting work. Nodal current kinetic energy uses annular weighted radial traces,
periodic axial averages, and explicit zero normal/azimuthal boundary velocities.
This is a declared interpolation whose physical convergence must be checked.

Integrated face watts combine full E cross B/mu0, electron/ion enthalpy with
factors beta/-alpha, and -delta (mu |w|²/2) J/e. The drift energy flux follows
the exact sum of species kinetic energy fluxes at zero barycentric flow.
The same nodal kinetic scalar enters gradient and flux. Global boundary power
vanishes under the bounded trace conditions. Local energy consistency is not
assumed exact; the result is never used to reset component energies.
