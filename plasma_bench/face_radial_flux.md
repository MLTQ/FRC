# Fluxes using independent radial momentum faces

Axial and angular momenta remain conservative cell means (n vz,n r vtheta).
Their face fluxes use actual stored radial q=r n vr at radial faces, not a
second interpolation of derived radial cell means. At axial faces radial
momentum is sampled from its same Qq physical lift. Density, both physical
pressures, angular geometry and enthalpy use the existing moment/nodal maps.

The inherited physical face integrals and directional divergence telescope
axial/angular momentum and the independent energy account. Radial momentum is
evolved by a separate variational equation and is not a third cell state here.
`radial_specific_momentum` exposes the exact transported traces for mechanical
work, avoiding a different reconstruction in the energy pairing. Impermeable
axis/wall and periodic z remain explicit. No exterior matter or diffusion.
