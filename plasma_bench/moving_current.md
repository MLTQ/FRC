# Moving-density finite-mass induction rate

Uses B*=B+(mu/e)curl(w), w=J/(ne), and actual supplied n_t and regular bulk
velocity. The full reduced equation is

    E = (delta w-v) cross B* - (mu/e) w cross curl(v)
        - (mu/e) w n_t/n + I_local J_t + Epressure_reduced
        + grad[phi_pressure + (mu/e)(v dot w-delta |w|²/2)].

This is equivalent to the original convective drift equation, including
w_t=J_t/(ne)-w n_t/n. The density and vorticity terms are not optional. The
instantaneous positive magnetic-space solver treats J_t implicitly. A complete
moving integrator must supply n_t from its own conservative continuity flux,
and evolve bulk momentum with physical drift stress and total energy with the
full species flux. This component alone advances neither n nor momentum.

Velocity has complete dual-node regular components (vr/r,vtheta/r,vz), with
zero radial wall velocity. Density/current/pressure retain the existing positive
CurrentThermal model and TraceHallCurl boundary contracts. Magnetic flux
inputs are validated for regular axis, fixed wall and periodicity before any
linear solve; these are bounded reference conditions, not plasma-vacuum data.
The returned full electric field restores both scalar potentials with the
commuting gradient. Nodal kinetic and
mixed potentials use the stated fixed-component current traces; this explicit
interpolation has a separately measured spatial error. No viscosity/diffusion,
density/current floor or state repair is introduced.

terms exposes the individual weak forcings for independent work and omission
controls. rate accepts an explicit optional electric source only for declared
forcing/verification; it is absent in unforced use. Returned fields are rates,
not a claim of a stable, conservative moving time step or vacuum coupling.
