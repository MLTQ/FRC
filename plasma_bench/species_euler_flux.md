# Species Euler flux in fixed velocity coordinates

This is a physical conservative face flux for each retained atomic level and
electrons, with all three momentum components and monatomic pressure. It does
not carry electromagnetic forces or enforce quasineutrality. The HLLC contact
construction follows the standard three-wave Euler solver; see the
[Clawpack Euler solver documentation](https://www.clawpack.org/gallery/pyclaw/gallery/shocktube.html).
The fallback is the conservative two-wave HLL flux, not a state correction.

Use the atomic moment scales: v0=sqrt(eV/m_neutral), p0=m_neutral*v0,
rho'=f*m_species/m_neutral, E in eV per reference nucleus, and p=(gamma-1)U.
Velocities w are relative to a fixed U, but the spatial cells remain fixed in
the laboratory. Normal laboratory velocity is u=w_axis+U_axis. Fluxes are

    F_f = f u
    F_P' = P' u + p e_axis
    F_E' = E' u + p w_axis.

The last formula includes the fixed-frame energy transformation; using
(E'+p)u would add spurious bulk-frame pressure work. Multiply all returned
dimensionless fluxes by v0 to obtain fluxes in the atomic state units times m/s.
Binding energy moves with the actual particle flux; the local produced-photon
and prescribed-heat ledgers are not material species and are not advected here.

Wave bounds include both endpoint acoustic speeds and positive Roe-averaged
thermal/velocity-variance contributions, without subtracting large bulk kinetic
energies. A vacuum side also includes its expansion-fan front speed. Both empty
sides have exactly zero flux and no atmosphere. HLLC star states follow mass,
momentum and energy jump conditions in the fixed frame. A vacuum face or negative
star thermal energy uses HLL instead. This adds explicit numerical diffusion;
it is not viscosity, conduction, a density floor or physical plasma edge physics.

Face number densities are scaled before contact-speed division, avoiding a
kilogram-mass product with a tiny population. Invalid actual input states or
nonfinite fluxes reject. General shock/vacuum-interface accuracy and plasma
boundary validity are not established by the smooth transport/decay examples.
