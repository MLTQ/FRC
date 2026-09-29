# Integrating the physics required by the target regime

The default Ar+ example has di/R about 4.10. The conditional local-edge scan
also crosses electron-inertia, nonlocal and carrier-drift scales. A single-fluid
scalar-eta run at those parameters cannot be presented as physically justified
formation. Choosing a denser/easier case does not discharge the original
low-pressure sandbox requirement; any alternate regime must be labeled.

## Next implementation and verification sequence

1. Add the missing axisymmetric toroidal magnetic and angular-momentum degrees
   of freedom, and derive the current/force/induction work pairing from parent
   species equations. Axisymmetry does not make the Hall term vanish in general.
   Initially verify a controlled fixed-domain limit before changing the moving
   FRC stress cases. Independent parallel-propagating circularly polarized
   Hall/ion-cyclotron and electron-inertial wave dispersion must check both
   handedness, phase speed, damping and energy, including time/mesh refinement.
2. Use separate electron and ion thermal energies, physical collision exchange
   and consistent current drift energy. The present current component is
   fixed-density/fixed-geometry: a moving extension additionally changes its
   mass matrix and kinetic-energy coefficient. Those variations contribute
   forces and transport and cannot be discarded while retaining its budget.
   Verify equal/opposite energy exchange and positive collisional heating
   independently; no unexplained total-energy repair is permitted.
3. State and check every reduction of the parent two-fluid momentum/continuity,
   species internal energy and Maxwell equations. Full parent momentum has
   `ms ns D_s vs = qs ns (E + vs cross B) - grad ps - div Pi_s + R_s`.
   Omitting displacement current, charge separation, pressure anisotropy or
   kinetic effects requires a measured local ordering. True vacuum has fields
   but no species mass/current; zero carrier density cannot be handled through
   a numerical gas floor or an arbitrary eta ramp.
4. Resolve which edge regime is actually supported: a collisional interior may
   require a nonlocal/kinetic or sheath interface before vacuum. Document the
   coupling conditions and energy/particle transfers, then construct an
   independent heated expansion reference. The original incompatible sharp-edge
   initialized FRC remains a stress case, not a secretly replaced target.
5. Couple prescribed spatial coils using MOVING_COIL_PLAN's geometric force and
   source-work reference. Establish chamber/end conditions and mesh/shock
   robustness, then run formation from unreversed plasma with a reproducible
   experimental interface and independent topology/force/validity criteria.

Step 1 now has a verified linear Hall/finite-mass reference (49 cylinder cases,
16 parallel modes) and a nonlinear compressible toroidal MHD extension (33 cases,
conservative angular momentum, full stress/energy and scalar resistive work).
See LINEAR_HALL_VERIFICATION and TOROIDAL_MHD_VERIFICATION. Nonlinear Hall is
now coupled and verified in COMPATIBLE_HALL_VERIFICATION (63 complete cases);
moving finite-mass and separate species thermal coupling remain missing.
TWO_TEMPERATURE_PLAN is the next bounded task: construct compatible charge flux,
electron pressure work and electron enthalpy transport, then verify warm species
dynamics and collisions. The current Hall solver has explicitly cold electrons.
The older incompressible oracle cannot be treated as an exact compressible mode.
The isolated dual electron-pressure/entropy component now supplies conservative
current flux, nonlinear Biermann/thermal exchange and discrete physical entropy,
with 46 complete measured records (ELECTRON_PRESSURE_VERIFICATION). It does not
yet advance moving ions or share HallMHD's boundary. Step 2 must next join primal
fluid and dual electron particle/energy fluxes without energy resets, then add
separate-temperature collision exchange and remaining physical transport.
Steps 2–5 remain required; neither wave nor nonlinear reference transport tests
validate the heated plasma boundary or actual formation.

The earlier 42-case Hall spatial operator failed its axis consistency gate
(21.46% at nr128), despite exact assembled work cancellation. It remains as a
historical control. The compatible successor corrects electric/energy spaces
and reaches 0.00830% axis error, then evolves all three compressible branches
and finite-amplitude variable-density flow with nonlinear implicit induction.
HALL_WORK_VERIFICATION retains the failure; COMPATIBLE_HALL_VERIFICATION records
the correction and actual time evolution. Neither reduces the required thermal,
finite-mass, heated-edge, vacuum or formation scope.

The already implemented linear `MaterialCurrent` limit remains useful for
coefficient, vacuum, source and electron-energy checks. Its 35-case records,
known coarse false reversal and thin-layer accuracy limits are retained in
CURRENT_VERIFICATION. It must not be installed as a complete moving two-fluid
solver merely because the component energy residual is small.

Parent equations and reduction caveats are discussed in
[NRL, species transport equations](https://www.govinfo.gov/content/pkg/GOVPUB-D210-PURL-LPS118567/pdf/GOVPUB-D210-PURL-LPS118567.pdf),
printed pp. 36–39, and
[Martin, generalized Ohm and current-energy derivation](https://ecommons.cornell.edu/server/api/core/bitstreams/17448f20-41b4-4457-a472-3ff41ab87529/content),
sections 2.2–2.3. These sources guide the derivation; their different experimental
geometries and approximations are not validation of this project.
