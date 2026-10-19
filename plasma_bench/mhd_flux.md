# Fluid fluxes and fast-wave bounds
`physical_flux` evaluates conservative poloidal ideal-MHD fluxes in either
direction. `rusanov` sets a common staggered normal B for both interface states
and recomputes face total energy at unchanged thermal pressure. Its local
Lax–Friedrichs jump term is explicit numerical diffusion, not physical viscosity
or resistivity. `fast_speed` supplies the exact ideal-MHD fast characteristic
bound, clamping only roundoff in its nonnegative discriminant. `ghost_cells`
reflects vr/Br at the axis and conducting radial slip wall, and wraps z.
`rusanov(..., magnetic_dissipation=True, magnetic_energy=True)` controls the
magnetic part of the energy flux:
- `magnetic_dissipation=False` removes the tangential magnetic-energy jump from
  the energy dissipation.
- `magnetic_energy=False` also removes the magnetic part of the central energy
  flux.

In either case the caller must supply the omitted transport through the CT
electric field's Poynting flux.
