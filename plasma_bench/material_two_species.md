# Material two-species free-boundary fluid

New moving fluid component using only MovingDomain connectivity/quadrature.
It does not call or replace the old P0 MovingPlasma solver. Material particles
at interior quadrature carry fixed particle weights and two separate specific
entropies. Density is particles/current cylindrical volume. An initial compact
profile can vanish exactly on the moving surface while all interior evaluations
remain positive. There are no exterior fluid cells, density/temperature floors,
pressure repairs or suppressed heating terms.

The canonical state contains general meridional node positions/momenta,
axisymmetric angular momentum and two material entropy arrays. The full physical
mass is mi+me. The consistent meridional mass matrix is fixed; the rotational
mass integrates M*n0*r² and changes with actual geometry. Its exact geometric
force is the centrifugal term at fixed angular momentum. Axis radial coordinates
and momenta are inactive; no outer position or normal velocity is prescribed.
The mass matrix is not lumped and velocities are not constrained to affine motion.

T_s[eV]=exp(s_s)*n^(gamma_s-1), U=sum_particles sum_s e*T_s/(gamma_s-1).
Pressure force is -dU/dx=integral p div(delta x)dV, using actual cylindrical
Jacobians. RK4 advances the full canonical equations; all stage/endpoint
Jacobians and positive interior species values validate without repair.
Kinetic energy, both internal energies and angular momentum are observable.
No numerical diffusion is introduced. There is no mesh remapping or shock model.

This component currently has no electromagnetic/current or collisional source
coupling. Its zero-current, equal-species-temperature expansion is an explicit
continuum fixture where Hall/current/Joule/exchange terms vanish, not a substitute
for the full finite-coil solver. General current-energy geometry forces, pressure
work, species exchange and vacuum/coil coupling must be connected before claiming
a complete moving plasma solver. Cold dilute physical-vacuum profiles are
mathematical fluid data, not a validated ionization/transport model for a lab edge.

Reference basis gradients are stored separately from the time-derivative method.
The initial direct-run attempt failed before stepping on a name collision; its
frozen source/failure remains in material_vacuum_sphere_r4_t800.

The material angle is now also evolved with theta_dot=Omega. Although the fluid
energy is independent of theta, a later magnetic pullback depends on its spatial
gradient. deformation(state) returns the full axisymmetric3-D material gradient
in orthonormal cylindrical frames, including r*grad(theta) and r/R. It enables
future full-vector field transport without silently discarding toroidal winding.
It does not itself implement that electromagnetic coupling. Earlier five runs
used identical fluid dynamics without storing angle; frozen sources are retained.

Six full trajectories and measured limits are in
docs/MATERIAL_CARRIER_EDGE_RESULTS.md. The zero-current expansion/rotation
component is exercised; complete moving electromagnetic coupling remains open.

MaterialIdealPlasma now composes this fluid with MaterialMagnetic forces/torques.
The fluid class itself remains unchanged; its original zero-current runs remain
valid. The new wrapper is only an ideal-MHD limit, not the missing full finite-mass
current/species/coil coupling. See MATERIAL_MAGNETIC_COUPLING.md for actual runs.

MaterialThermalCurrent now composes the unchanged fluid with finite-mass current,
general species entropy/pressure work, prescribed resistance/heat and conditional
exchange in the meridional/toroidal-field sector. Read MATERIAL_THERMAL_COUPLING.md
for complete coarse runs, the preserved fine entropy failure and ongoing upwind
replacement. Poloidal/vacuum/coil and physical edge validity remain incomplete.
