# Spatial transport of reacting species moments

`AtomicTransport` gives the actual atomic moment state conservative axial fluxes
for number, all three momenta and total translational energy of every heavy level
and electrons. Cell values are densities divided by ONE fixed reference nucleus
density Nref. Nucleus density can vary and advect; cells are never renormalized
to unit nucleus count. The existing source uses Nref*f_e as electron density and
accepts these extensive normalized numbers without changing its reaction laws.

The cells are fixed in the laboratory and currently periodic. Fixed velocity
coordinates follow atomic_frame.md; the corresponding pressure work is included
in species_euler_flux.md. No particle velocity or mass is discarded, and no
charge/current projection is made. The flux transports bound internal energy
through the actual populations. Produced photons and prescribed input heat are
local bookkeeping densities with their source rates; photon propagation is not
implemented by pretending they move with matter.

Spatial reconstruction uses componentwise minmod slopes of number, pressure and
relative velocity. HLLC resolves the material contact; HLL handles vacuum/star
failures. Both introduce numerical diffusion, strongest near limited extrema and
discontinuities. This is explicitly numerical, not physical heat/neutral transport.
SSP-RK2 advances fluxes; source-half/transport-full/source-half composes the local
reaction or reacting/colliding model. Source heat samples the true quarter and
three-quarter times. There is no energy repair, population floor or atmosphere.
Every stage and endpoint validates the original actual species thermal state.

`timestep` gives the explicit acoustic CFL step, including the physical electron
mass. Reaction/collision accuracy can require an additional user-specified dt.
`totals` measures periodic integral nucleus, mass, charge, vector momentum and
frame-energy + binding + produced photons - input heat. Fixed-frame invariants
imply laboratory invariants; no total is imposed on the evolved state.

This is a spatial operator component, not a replacement for the moving-coil
solver. No Lorentz/electric constraint, radial geometry, sheath, open lab ends,
photon propagation or free-boundary interface tracking is present. Numerical
charged transport can separate charge unless electromagnetic feedback enforces
the appropriate equations. The exact contact fixture has zero current/charge
and constant pressure; the decay fixture has exactly zero charged particles.
Neither is a physical FRC formation experiment. Those missing couplings remain
required before applying arbitrary charged spatial evolution to the laboratory.

Completed neutral expansion shows nonzero numerical density beyond the exact
vacuum front on every measured grid, despite zero initial exterior matter.
This operator must not replace the formation solver's tracked free boundary.
Measured histories and convergence/limitations are in docs/ATOMIC_TRANSPORT.md.
