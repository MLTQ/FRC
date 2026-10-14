# Moving two-species fluid at a vacuum edge (2026-10-14)

This is progress on frc-vacuum.3, not a completed moving electromagnetic plasma
solver. The existing finite-coil solver remains available and its evidence is
unchanged. The new MaterialTwoSpecies component releases the fixed fluid
boundary while carrying actual physical mass, two entropies and angular momentum.
The completed trajectories below are the zero-current, equal-temperature Euler
limit; they do not validate a heated, current-carrying plasma edge.

## Physical change

Material quadrature carries particle weights dN and two specific entropies.
At each stage, n=dN/dV uses the actual moving cylindrical volume. Both species
obey T_s=exp(s_s)*n^(gamma_s-1), with gamma_s=5/3 and full mass mi+me.
The meridional mass matrix integrates (mi+me) Ni Nj dN. Angular inertia integrates
(mi+me) r² Ni Nj dN and changes with position. Canonical angular momentum gives
the associated centrifugal force; the outer velocity/shape is never prescribed.
Pressure force is the exact discrete derivative of both species internal energy.

The body has general Q1 meridional node motion, not an affine restriction.
Only the symmetry-axis radial coordinate is constrained. This reuses the old
MovingDomain connectivity and quadrature, not its P0 MovingPlasma dynamics.
There are no exterior fluid cells, density/temperature floors, numerical
viscosity, energy correction or silently suppressed heat/current terms.
The compact reference carrier profile supplies a vanishing density trace;
interior quadrature samples stay positive. General arbitrary input arrays are
not a guarantee of that boundary trace: the physical reference profile must be
specified. No field is reconstructed by extrapolating a positive density floor.

Material angle is now a stored coordinate with theta_dot=Omega. The fluid energy
is independent of angle, but a later full magnetic pullback depends on its
spatial gradient. A full axisymmetric deformation-gradient helper includes
r*grad(theta) and r/R. Electromagnetic transport using that helper is not yet
implemented or validated. Pressure/current/current-energy geometry forces and
species relative transport must still be connected.

## Independent continuum baseline

The reference belongs to the compressible affine physical-vacuum Euler family
in [Sideris, Global Existence and Asymptotic Behavior of Affine Motion of 3D
Ideal Fluids Surrounded by Vacuum](https://web.math.ucsb.edu/~sideris/pdffiles/explicit-12-09-15.pdf).
This is a mathematical fluid reference, not a plasma-transport, stability or
confinement result. The normalization below follows directly from this run's
physical two-species EOS.

For initial ellipsoid axes a0,c0 and reference spherical coordinate s:

    n0=nc*(1-s²)^(3/2), Ts0=Tsc*(1-s²), s<1; no matter outside.
    C=5e*(Tec+Tic)/(mi+me), D=a²c/(a0²c0).
    a''=C*D^(-2/3)/a+Omega0²*a0^4/a³,
    c''=C*D^(-2/3)/c, Omega=Omega0*(a0/a)².

For a sphere at initial rest, R(t)=sqrt(R0²+C*t²) is an exact closed form.
The elongated rotating reference uses an independently integrated small ODE,
not the production node-force equations. Total particles N=pi²*a0²*c0*nc/8,
U=N*e*(Tec+Tic)*15/16*D^(-2/3), and
K=(mi+me)*N*(2*a_dot²+c_dot²+2*a²*Omega²)/16.

Runs use a0=35mm,c0=35mm (sphere) or70mm (rotating), nc=1e20/m³,
Te=Ti=1eV, Ar+ physical mass, and800 steps through t=a0/sqrt(C), about7.12us.
Rotating data use Omega0=0.5*sqrt(C)/a0. B=J=0 and equal temperatures make
Hall/current stress/Joule heat/Coulomb exchange vanish in this fixture. A cold
edge is not a justified fully ionized laboratory transport/ionization closure.

## Completed trajectories

Three spherical and two rotating runs completed. Actual batch90694 and87440
exits0, including all child exits, were collected. The initial direct r4 attempt
failed before stepping on a method/array name collision; its source and failure
are preserved in material_vacuum_sphere_r4_t800. Corrected r4 has suffix_v2.
All accepted states, canonical momenta, entropies and frozen source closures are
saved in the corresponding outputs/material_vacuum_* directories.

Final errors are normalized by initial radial axis, speed sqrt(C), or the
current central density/temperature as indicated:

| Mode/grid | Maximum surface error / a0 | Velocity RMS / sqrt(C) | Physical-space density RMS / center | Temperature RMS / center |
|---|---:|---:|---:|---:|
|Sphere4/8|0.0197782|0.00558611|0.0133015|0.0147663|
|Sphere8/16|0.00477359|0.00139297|0.00332579|0.00371430|
|Sphere16/32|0.00114279|0.000347990|0.000831438|0.000929990|
|Rotating8/16|0.00282588|0.00123271|0.00325642|0.00390218|
|Rotating16/32|0.000698788|0.000307962|0.000813368|0.000977048|

Refinement gives roughly second-order physical convergence. The finest final
surface errors are40.0um for the sphere and24.5um for the rotating ellipsoid.
Angular-speed RMS error / (sqrt(C)/a0) falls0.000648514→0.000162403.
Canonical total angular momentum is constant by the rotational symmetry; the
independent changing angular-speed comparison tests its physical evolution.

Native energy change remains below4.0e-14 relative throughout these trajectories.
That conservation does not remove geometric accuracy error: the polygonal
initial mass is3.8067%,0.96086%,0.24078% low on the three meshes. Initial
internal/kinetic energy discrepancies are also separately recorded. Physical-
space density/temperature comparisons include the initial geometric mismatch;
material-coordinate errors are recorded separately. No energy was repaired.
Peak RSS for the five successful runs stays below123MB; elapsed seconds in table
order are5.688,20.224,78.699,22.054,85.227. No unit tests or existing suite ran.

The additional rotating16/32 trajectory storing material angle also completed
all800 steps. Actual session12290 exited0; elapsed84.966s,peakRSS127369216B.
Its final phase RMS error is0.000100835rad versus reference phase0.374825535rad.
The physical boundary/velocity/energy equations are unchanged from the earlier
run; its independently stored angle is required for a fully material rotational
map. The new deformation helper has not yet been exercised by a coupled field
trajectory. All six successful runs and the initial setup failure are terminal;
no scientific process remains live.

The inspected figure is
outputs/material_vacuum_comparison/material_vacuum_expansion.png (with PDF from
the same figure). Observer95156 completed with actual exit0. It shows actual
surfaces, spatial error and decreasing rotation rate from saved states.

## Next physical integration

Keep frc-vacuum.3 and the complete FRC goal open. The next work is full moving
magnetic/current coupling on this material state, including toroidal winding,
current kinetic-energy geometry derivatives, finite-coil mutual/source work and
both species' relative transport. A carrier-weighted current-velocity space is
a candidate way to make J vanish with n and avoid a hidden inverse-density
floor; it is not yet an adopted or verified discretization. It must retain the
original physical inertia, Hall, advection/stress and pressure work equations.

Read SURFACE_COMPATIBILITY.md before enabling edge heating. Passing these
zero-current expansions does not resolve a current-carrying heated boundary,
kinetic/sheath or ionization effects, or FRC formation. Do not replace the full
finite-coil model with this fluid component or the older scalar moving solver.
Continue with actual coupled physical trajectories, not another test framework.
