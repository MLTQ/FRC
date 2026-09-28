# Fluid fluxes and fast-wave bounds
`physical_flux` evaluates conservative poloidal ideal-MHD fluxes in either
direction. `rusanov` sets a common staggered normal B for both interface states
and recomputes face total energy at unchanged thermal pressure. Its local
Lax–Friedrichs jump term is explicit numerical diffusion, not physical viscosity
or resistivity. `fast_speed` supplies the exact ideal-MHD fast characteristic
bound, clamping only roundoff in its nonnegative discriminant. `ghost_cells`
reflects vr/Br at the axis and conducting radial slip wall, and wraps z.
