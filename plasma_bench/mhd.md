# Developmental dynamic MHD solver
`AxisymmetricMHD` evolves conservative mass, both poloidal momenta, total energy
and corner flux together with SSP RK2. `MHDParameters` exposes gamma, physical
scalar resistivity and CFL; no diffusion/stability term is hidden. Timestep
restrictions combine directional fast waves and explicit resistive diffusion.
`step` validates both stages before committing; no pressure/density floor or
energy repair is permitted. `advance` retains partial state if a step budget
fails. `diagnostics` reports conservation, energy partition, physical heating
estimate and divergence. Fluid fills the domain; z is periodic and the radial
boundary is a slip conducting wall. This is not yet the chamber/vacuum formation
solver. The original experimental Stage 2 gate remains separate until coupling
is verified. See docs/MHD_NUMERICS.md for equations and numerical dissipation.
