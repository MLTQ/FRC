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

**2.5D face states (2026-10-20).** Eight-column states (…, Ω = v_θ/r,
g = B_θ/r) passed with a face `radius` use `toroidal_flux.physical_flux`:
conservative angular-momentum flux r(ρ v_θ v_n − B_θ B_n/μ0) and B_θ/r flux
g v_n − Ω B_n, with the fast speed including B_θ. Both magnetic-energy switches
then remove the full B² (including B_θ²) and v·B (including v_θ B_θ).
Six-column states follow the unchanged poloidal path.
