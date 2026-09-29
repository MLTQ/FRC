# Moving axisymmetric parent-species manufactured reference

Defines analytic positive n(r,z,t), three velocity components, both nonuniform
species pressures, poloidal/toroidal B, its exact curl current and a Faraday-
consistent E. Explicit particle fluxes determine n_t=-div(n v), with zero
radial flux at the cylinder wall and periodic z. Because div J=0, both species
have the same continuity law. Regular radial/azimuthal fields vanish at the
axis. Psi stays fixed at axis/wall and Btheta(R)=0.

The velocities vi/ve and parent species kinetic/enthalpy flux are assembled
directly. Complex-step derivatives of analytic fields yield the actual force
and thermal sources in each species momentum/internal-energy equation. They
are explicit manufactured forcing, not an unforced equilibrium or a formation
experiment. Finite differences independently audit these derivatives.

Different gamma_e/gamma_i, changing density, bulk rotation/compression/shear,
non-Beltrami current and finite electron inertia exercise the full reduction.
The default low-density fixture exposes inertia while retaining physical Ar+
and electron masses. The original n=1e20 regime is a separate control. Local
cylindrical differential operators evaluate at r>0; axis limits are checked
by approaching the axis and by regular production interpolation tests.
