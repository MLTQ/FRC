# Architecture and stack decision

Python 3.11+, float64 NumPy arrays, SciPy elliptic integrals and sparse LU,
Matplotlib/ContourPy visualization and topology, pytest verification.
The initial environment is macOS ARM64 / Python 3.13, not the eventual RTX 4090
machine. CPU correctness comes first; no accelerator is assumed or required.

NumPy/SciPy provide the required special functions and a sparse direct solver
without writing numerical infrastructure. Geometry-dependent loop bases and
diffusion factorizations are reused. For this modest axisymmetric problem, a GPU
transfer/JIT layer would add complexity before we have a performance bottleneck.
JAX can support future differentiable/JIT array kernels, Numba can accelerate
profiled local loops, and PETSc becomes attractive for larger implicit solves.
None is required by a public interface. Revisit after measured profiling and the
Stage 2 solver choice, on the actual target GPU operating system.

The bounded [external-engine comparison](ENGINE_COMPARISON.md) now verifies an
optional WarpX 26.09 CUDA installation on the remote RTX 4090. It is isolated
from the core dependencies and is not a production backend. Experiment/source
definitions and result diagnostics remain reusable; adopting kinetic dynamics
still requires explicit model, vacuum/edge and output-conversion contracts.
The current custom axisymmetric research path is retained.

`current_thermal`, `current_thermal_midpoint` and `current_thermal_flux` form a
separate constrained fixed-n, zero-barycentric-flow subsystem with physical
reduced-mass current inertia, both species entropy and the full drift-energy
flux. Stable current increments preserve weak signals under a large background;
stored-state backward error is distinct from nonlinear Ohm error. No moving
solver was silently switched to this component. Its analytic references,
baseline archives, independent energy/current audit and precision controls are
documented in CURRENT_THERMAL_NUMERICS and CURRENT_THERMAL_VERIFICATION.

References: [SciPy cached factorization](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.factorized.html),
[JAX supported installation platforms](https://docs.jax.dev/en/latest/installation.html).

```
geometry + typed configuration
      |         waveforms -> coil definitions
      |                          |
      +---- analytic vacuum field basis
      |                          |
      +---- plasma solver protocol (Stage 0 / Stage 1 / future MHD)
                                 |
                     result arrays + metadata
                                 |
                    diagnostics / topology
                                 |
             plots / time slider / sweeps / saved NPZ
```

Modules have one responsibility and adjacent .md contracts. Arrays use `(nr,nz)`
and histories `(nt,nr,nz)`. SI fields and axis conventions are shared; coils do
not import plasma solvers. Frozen dataclasses validate geometry, material and
source configuration. Waveform objects supply value and derivative; callbacks
are script-only and rejected by reproducible serialization.

Results include full configuration, source/version fingerprint, dependency
versions, actual time coordinates, warnings, units, boundary conditions, model
stage and numerical method. Saved archives are numeric NPZ without pickle;
metadata is JSON. JSON configurations and a CLI keep runs reproducible without
a notebook. Sweeps replace typed configurations, not hidden globals.

The developmental dynamic kernel is `mhd.py`, with explicit modules for annular
mesh geometry, conservative state, Rusanov fluxes, compatible induction,
axis-regular reconstruction and analytic benchmarks. It evolves density,
poloidal momentum, total energy and corner flux with SSP RK2. Both first-order
reference and default second-order spatial errors are measured. The independent
`mhd_equilibrium` module supplies a continuum closed-flux Grad–Shafranov fixture;
it does not use the fluid solver to construct its reference force balance.
The `mhd_positivity` helpers bound density/thermal reconstruction traces while
leaving conserved cell averages intact. The solver records accepted-stage
activation counts; trace admissibility does not imply a full-update theorem.
It fills a verification cylinder with fluid, with a conducting radial wall and
periodic z. The original Stage 2 formation entry point still refuses execution
until moving-plasma/vacuum and source-work coupling are verified; it does not
silently substitute these verification boundaries for the physical experiment.

The separate `coupled.py` reference integrates moving gas in a fixed chamber
with prescribed analytic coils and a magnetic-only exterior. `coupling_mesh`
maps fluid cell velocities/heat to current-bearing field nodes; `coupled_hydro`
advances conservative gas transport under reflecting walls/ends; and
`coupled_magnetic` pairs force/advection by adjoints in an implicit midpoint
solve with exact elliptic vacuum rows. Strang splitting combines the steps.
All stages must succeed before state, time or cumulative work is committed.

`CoupledState.gas` stores internal plus kinetic energy, while the global field
owns induced self-energy. This deliberately differs from `MHDState`'s total
energy; converting arrays by simple copying would double-count magnetic energy.
Diagnostics name the response-energy budget explicitly and retain source work,
physical heat, vacuum current and nonlinear residuals. `coupled_benchmarks` and
`vacuum_benchmarks` supply independent mathematical references; separate example
measurement/report modules save complete configurations, histories and source
fingerprints. See `COUPLED_NUMERICS.md` and `COUPLED_VERIFICATION.md`.
The fixed chamber reference does not yet supply the moving free plasma interface,
physical transport or boundary validity required by the formation workflow.

`open_boundary` adds an analytic decaying multipole exterior to the finite field
grid. Its modal amplitudes follow interior nodal unknowns in the global energy
matrix, with algebraic vacuum stationarity rows. No fluid or resistivity is
assigned to modes. Boundary nodal values are their trace; the entire exterior
magnetic energy is counted. `CouplingMesh.interior/full` pack and expand induced
unknowns, while `grid_vector` handles arbitrary physical grid data such as bias
and coil flux without projecting them into the decaying basis. Diagnostics keep
modal stationarity separate from actual vacuum current density. Independent
asymmetric loop-integral and dynamic pulse checks are in `open_benchmarks` and
the open-boundary report. Fixed-boundary comparisons remain selectable.

`conservative_coupled` is the unsplit conservative/CT comparison using the same
magnetic vacuum. `vacuum_flux` builds shared fluid faces and corner electric
fields; `vacuum_trace` supplies exterior traction and explicit vanishing
numerical interface diffusion. `flux_update` owns geometric flux assembly and
shared local replacement. `flux_limiting` detects complete Euler admissibility,
expands only the affected region, and retains a global comparison mode. Exact
fluid energy-flux conservation and independently measured whole-space response
error are separate contracts; neither pressure floors nor energy resets hide
the latter. See the conservative coupling and local-limiter reports.

The separate experimental `moving_*` modules evolve a material surface rather
than restricting fluid to a fixed chamber. `moving_geometry` owns material Q1
cells, massless exterior mesh motion, a curved matching sphere and integration.
`moving_field` solves stationary vacuum degrees of freedom and derives magnetic
geometric force from full-space induced energy. `moving_advection` supplies
fixed nodal flux labels and an explicit regular-axis stretch approximation.
`moving_plasma` owns actual material positions/velocities, fixed masses and
adiabatic cell entropy, with second-order Verlet and prescribed uniform-bias
work. No exterior fluid arrays are allocated. Dynamic field arrays use node
and element ordering, not the earlier rectangular-grid layout.

This reference is ideal and smooth-flow only: it does not silently replace the
resistive fixed-chamber solvers or unlock Stage 2. Its accepted state plus fixed
material constants and full magnetic coefficients are archived by independent
`examples/moving_*` measurement/analysis modules. Solver/measurement provenance
is frozen during a study; report/plot hashes are separate. The measured flux
interpolation defect, magnetic translation-force error and equilibrium drift
are part of the numerical contract, alongside missing shock/remapping and
transport support. See MOVING_PLASMA.md and MOVING_PLASMA_VERIFICATION.md.

`material_resistivity` adds a consistent plasma-only r³-weighted mass matrix,
weak Jtheta/r recovery and implicit midpoint diffusion against the same
full-space induced-field stiffness. Vacuum has no mass or resistivity. Physical
eta J² supplies nonnegative cell heat. `resistive_substeps` uses step doubling
to control magnetic, current and local-heat differences, accepting the two-half
solution without extrapolation or energy repair. Raw midpoint can falsely
reverse an underresolved decaying mode even with exact energy bookkeeping;
independent stiff controls retain that failure and measure adaptive improvement.

`resistive_moving.ResistiveMovingPlasma` composes resistive half-steps with the
unchanged ideal Verlet step. It evolves flux labels and cell entropy, preserves
fixed mass, and commits all stages atomically. Eta=0 retains the ideal path.
Current, cumulative cell heat, labels and entropy join the archived geometry,
velocity and fields. Physical scalar eta is a prescribed input, not a transport
closure. Uniform bias is still the only moving-source geometry. The completed
13-case continuum and 20-case moving studies distinguish current/heat convergence
from the large remaining flow error; see RESISTIVE_MOVING_NUMERICS.md and
RESISTIVE_MOVING_VERIFICATION.md. Neither this extension nor its local tolerance
unlocks Stage 2 or replaces mesh, surface thermodynamics and regime validation.

Optional radial grading in `MovingDomain` resolves material and magnetic-vacuum
shells near their shared surface without changing physics; the default uniform
coordinates are retained exactly. Angular resolution remains independent.
`examples/surface_sampling` inverts the initial geometry to compare common
physical material points. `material_overlap` intersects initial cells and uses
positive quadrature, preserving discontinuities in the comparison rather than
smoothing thermal fields. Cached paired-order integrations retain exact input
and analysis hashes. Independent moments, gas expansion and sphere decay verify
the added geometry/control, separately from actual surface-flow convergence.
See SURFACE_NUMERICS.md and the conditional thermal trace audit in
SURFACE_COMPATIBILITY.md; mesh refinement cannot by itself validate missing
physical interface assumptions.

## Current inertia and conditional transport

`current_inertia.py` is a separate fixed-material linear current component. It
reuses the existing full-space stiffness/plasma mass, adds positive carrier
inertia and relaxation matrices, and solves a mixed Ampere/Ohm increment system.
Magnetic energy, electron drift energy, per-cell heat and source work remain
separate. No incomplete current-energy force or Hall term is inserted into
`ResistiveMovingPlasma`. Three independent continuum references, time order and
physical-volume energy checks are recorded in CURRENT_VERIFICATION.

`coulomb_transport.py` supplies conditional classical coefficients and local
physical-scale diagnostics. `diagnostics.characteristic_scales` includes these
as a named initial reference without changing conductivity. Zero-density
vacuum has no transport state; invalid Coulomb closure is not floored. The
next Hall/toroidal and electron/ion thermal integration contract is
EXTENDED_PLASMA_PLAN. Stage 2 remains gated on full formation requirements.


## Linear axisymmetric Hall reference

`hall_parameters` defines immutable SI species coefficients and explicitly named
massless comparison models. `hall_modes` owns the paired four-component linear
magnetic/velocity dynamics with positive current-energy mass and implicit
midpoint. `hall_cylinder` supplies numerical radial FE matrices, an odd periodic
Fourier grid, regular poloidal/toroidal reconstruction and physical-volume
quadrature. It evolves toroidal field and angular momentum in a bounded linear
incompressible limit; these arrays are not compressible `MHDState` arrays.

Independent Bessel/circular references, a separately constructed spatial
exponential and reproducible measurement/report modules record continuum and
time errors separately. The field archive includes every initial/final mode,
not only plots. All 49 cylindrical and 16 parallel modes are complete. No
vacuum boundary or thermal closure is inferred from these reference walls.
LINEAR_HALL_NUMERICS and LINEAR_HALL_VERIFICATION define the contract; the next
compressible conservative integration is TOROIDAL_INTEGRATION_PLAN.


## Compressible toroidal extension

`toroidal_state` defines an explicit six-variable state with conservative axial
angular momentum and Btheta/r, retaining exact regular annular r² moments for
rotation/magnetic energy. `toroidal_flux` and `toroidal_reconstruction` provide
full three-component ideal fluxes, speed bounds and centrifugal/hoop stress.
`toroidal_induction` pairs a compatible positive radial current mass with
resistive Faraday, Poynting flux and heat; it includes the finite axis electric
field. `toroidal_mhd` atomically advances the coupled nonlinear compressible
state with SSPRK2. The old four-variable MHD API is unchanged, and the new
zero-toroidal RHS agrees with it at both numerical orders.

Independent references include exact finite-amplitude translating rotation and
conducting Bessel decay. Physical-volume errors of the limited reconstruction
are separate from cell-average comparisons. The complete 33-case records retain
all conservative arrays, sources and numeric-column histories. The reference
cylinder still has periodic z and conducting slip radial wall; it has no Hall,
free-vacuum, neutral or separate species-thermal coupling. TOROIDAL_MHD_NUMERICS
and TOROIDAL_MHD_VERIFICATION define the measured scope. NONLINEAR_HALL_PLAN
specifies the next current/induction work pairing and finite-sound-speed oracle.

## Historical experimental nonlinear Hall spatial component

`hall_curl` differentiates the existing magnetic state energy and pairs its
adjoint Faraday map with a common electric quadrature. `hall_induction` projects
the full nonlinear J cross B/(ne e) with the same masses, exactly canceling
assembled Hall work. `hall_poynting` independently evaluates face energy
transport and records the local numerical mismatch without repair or heating.
There is no new fluid force, time integrator or connection to ToroidalMHD.

The finite-pressure reference derives all compressible branches independently;
the nonlinear analytic reference includes variable density, all field/current
components and regular-axis limits. Forty-two spatial cases show global
convergence and roundoff Hall work, but expose a persistent 21.46% first-annulus
gdot error. This explicitly prevents integration. HALL_WORK_NUMERICS and
HALL_WORK_VERIFICATION describe the experimental stencil and failed gate;
NONLINEAR_HALL_PLAN requires a consistent regular-axis energy/current/electric
revision before local transfer and actual stable nonlinear wave evolution.

## Compatible compressible Hall evolution

`compatible_magnetic` owns an explicit positive poloidal energy for bilinear
psi in x=r²/2,z; `hall_compatible` supplies matching electric/magnetic sampling
with a commuting physical Faraday curl. `hall_matrices` builds the sparse
magnetic Hessian, skew Hall and positive resistance forms. `compatible_poynting`
integrates face work with the same axis moments. `hall_midpoint` solves their
nonlinear midpoint equations with checked actual linear/nonlinear residuals.
Its LU cache is a preconditioner, not an approximation substituted for the
changing physical matrix. Failed trials leave physical input intact.

`hall_fluid.HallMHD` composes ideal-fluid half steps and full implicit magnetic
steps; all fluid momentum remains in the existing full Maxwell-stress flux.
`ToroidalState.magnetic_quadrature` explicitly carries the energy convention
through pressure, stages, copies and diagnostics; the legacy default stays
unchanged. No state energy conversion, numerical atmosphere or Hall heat exists.

Independent references, protocols, records, spatial/evolution measurements,
supplemental controls, physical-energy audit, comparisons and figures are split
among the `examples/compatible_hall_*` modules. All 54+9 cases have complete
source/configuration/archive/history checks. The report separates global work,
first-order physical field representation, actual mode damping/phase, local
transfer, nonlinear mesh differences and second-order time convergence.
COMPATIBLE_HALL_NUMERICS/VERIFICATION define this bounded reduction; the next
charge/enthalpy and electron/ion thermal contract is TWO_TEMPERATURE_PLAN.

## Electron pressure and thermodynamic transport

`dual_current` derives locally conservative integrated charge fluxes directly
from toroidal magnetic work, including the physical wall trace and Poynting
boundary term. `thermal_means` implements stable logarithmic means with
independent high-precision checks. `electron_pressure` uses the same entropy
mean in electron drift and pressure Ohm; `electron_midpoint` pairs magnetic
midpoint with a logarithmic temporal state to preserve discrete energy and
physical electron entropy. Failed trials preserve inputs; no state repair,
temperature floor, ad hoc wall return current or extra heat is introduced.

Electron scalars occupy full-domain nodal dual volumes. They are not silently
mapped onto the existing primal fluid energy. This fixed-density pressure split
has explicit constant toroidal-field/density wall traces, generally nonzero
local pressure Ez and zero integrated boundary power. That boundary differs
from the prior conducting Hall wall, so the component is not inserted into
HallMHD. The next required integration is conservative primal/dual mass,
electron pressure/energy and boundary work; TWO_TEMPERATURE_PLAN retains it.

`electron_reference`, `electron_cases`, `electron_records` and
`electron_measurements` retain independent physical oracles and all 46 records.
`electron_audit` replays actual physical pressure/energy/entropy, checks every
archive/history, and separates time and nested-state norms. The report plots
continuous pressure error and actual battery energy transfer. Equations and
limits are in ELECTRON_PRESSURE_NUMERICS and ELECTRON_PRESSURE_VERIFICATION.

## Moving ideal species temperatures

`dual_mass` maps primal ion density and its actual integrated mass flux to the
dual thermal grid. The overlap and flux divergence commute locally; the reverse
Q1 scalar average is its volume adjoint. `ToroidalMHD._rhs_with_fluxes` exposes
the shared mass flux while preserving its previous RHS contract and arithmetic.
`bulk_entropy` transports electron U with that flux; `two_temperature` recovers
pe from changing n and pi from conservative total thermal energy. `warm_mhd`
integrates both through the same two ideal RK stages with atomic validation.

The total fluid energy is never reset to fit the thermal partition. This ideal,
equal-gamma, common-velocity path reproduces the prior total fluid trajectory.
It rejects resistivity until heating partition is derived; it has no Hall,
electron-current enthalpy, pressure battery, collisions, conduction or finite
electron mass. Positive cell temperatures and exact energy do not guarantee
physical entropy monotonicity or accurate shock partition.

Independent contact and finite-amplitude simple-wave references, 30 main cases,
six stronger-gradient/common-time controls, independent energy replay and
continuous pressure/temperature quadrature live in separate `examples/warm_*`
modules. WARM_ADIABATIC_NUMERICS/VERIFICATION document errors and the remaining
entropy defect. The next gate connects electron-current enthalpy/Poynting work
and Hall with a compatible boundary, as specified in TWO_TEMPERATURE_PLAN.

## Exact finite-mass moving-current contract

`finite_mass_moments` implements the exact species reduction, barycentric mass,
drift stress, current energy and full bulk/current energy flux. The separate
fixed-density `current_thermal` path retains its previous source identity.
`dual_kinematics` supplies regular nodal bulk velocity and cylindrical vorticity;
`moving_current` adds moving-density and bulk-gradient terms to instantaneous
induction. `inertial_induction` eliminates the mixed equations into a positive
magnetic-space solve with checked original residuals and explicit linear
iterative refinement. No moving fluid time integrator is hidden in these modules.

Independent analytic parent-species fields and manufactured force/heat live in
`examples/moving_species_reference`; the continuum audit checks both species
and local total energy. Separate case, work, measurement, provenance, audit,
plot and report modules preserve nineteen full rate records and the initial
failed attempt. MOVING_CURRENT_NUMERICS/VERIFICATION specify equations, domains,
measured errors and the still-required discrete momentum/thermal integration.
Actual conservative continuity must supply n_t, and a consistent primal/dual
velocity and stress/energy flux must precede full moving finite-mass evolution.

## Conservative moving finite-mass integration

`finite_mass_state` now defines an explicit total-species-mass conservative state
and its magnetic/current/bulk/electron/ion energy partition. `primal_velocity`
maps mass-weighted regular bulk velocity onto dual nodes. `finite_mass_flux`
and `finite_mass_transport` separate Rusanov bulk transport from full computed
Poynting power, drift stress and mixed current/enthalpy transport. The shared
mass flux supplies the actual n_t used by `moving_current` and electron entropy.

`finite_mass_rhs` evaluates the entire derivative without fixed-n splitting;
`finite_mass_mhd` advances it by explicit RK4 with atomic stage validation and
source/boundary energy accounting. A positive filled periodic cylinder is the
only implemented domain for this integrator. There is no physical vacuum edge,
collisional/neutral transport, general stability theorem or shock entropy rule.

Independent manufactured parent species and exact finite-mass compressible
waves live in separate `examples/finite_mass_*` references. Measurement, physical
energy/entropy audit, time comparison, source provenance and plots remain separate
modules. The 21 baseline rates preserve their original snapshot; 18 trajectories
have an expanded measurement fingerprint with unchanged production sources.

`annular_polynomial` is a separately verified diagnostic for nearly solenoidal
mass-flux cancellation. The original FiniteMassMHD does not select it; the
separate shared smooth transport described below uses its interpolation helpers.
FINITE_MASS_NUMERICS/VERIFICATION retain pressure errors rather than treating
global conservation as proof of accurate temperatures or formation physics.

## Shared smooth finite-mass reference

`annular_cell_polynomial` supplies transverse polynomials whose integrals recover
the actual input averages. `smooth_finite_mass_flux` uses those with shared cubic
face values for mass, momentum, enthalpy and kinetic transport; `smooth_bulk_entropy`
advects conservative U using the ratio of separately reconstructed U and n with
the actual mapped particle flux. `smooth_finite_mass_rhs` assembles the full
unsplit derivative; `smooth_finite_mass_mhd` inherits the original RK4, positive
stage checks and atomic budgets. The limited baseline remains a separate source.

`kernel_induction` algebraically separates exact separable curl kernels and
keeps the strict small-signal solve, restoring full E for work. Original full
equation backward errors and descriptive rate-relative conditioning are distinct.
This module does not supply a general Helmholtz decomposition or cure the
slow-mode spatial cancellation defect in velocity/current representation.

Independent contact parent equations, projection, continuous physical norms,
slow-mode oracle controls, precision-failure records, measurement protocols,
source archives, read-only audits and plots remain separate `examples` modules.
All 31 rates/27 trajectories are preserved, including old failed serialization
evidence. SMOOTH_FINITE_MASS_NUMERICS/VERIFICATION document improved pressure
transport alongside unresolved magnetic/energy errors and entropy defects.
No shock scheme, physical transport closure, exact-vacuum formation gate or
3-D stability/confinement result is enabled by this smooth bounded reference.

## Higher-order compatible fields and affine harmonic storage

`compatible_polynomials` supplies local nodal interpolation and interval-average
edge functions whose derivatives commute with incidence. `compatible_geometry`
constructs regular cylindrical B/E sampling and physical electric/magnetic
Gram matrices from separable tensor factors. Its fixed-axis/wall mixed trace
retains the existing Faraday map; no general boundary or vacuum substitution is
implied. `compatible_velocity` reconstructs degree-3/5 physical velocities from
actual conservative annular averages, without fitting them to current.

`compatible_current` supplies the physical constant-density mass/inertia/Ampere
contract consumed by the unchanged instantaneous inertia solver. `harmonic_magnetic`
stores a constant axial bias and increment separately, preserving weak currents
without discarding full physical field/local cross energy. Its rounded-total
export is explicitly lossy and is not an evolution input.

Independent root-product field norms, all-branch wave protocols, bias precision
diagnosis, separate source provenance, original-equation audits and scientific
plots live in `examples/compatible_*`. Both the first rounded failure and final
affine study are preserved. COMPATIBLE_FIELD_NUMERICS/VERIFICATION define what
has been proved and the remaining common work/particle-transport contract.
The nonlinear SmoothFiniteMassMHD remains unchanged: variable-density inertia,
current energy, thermal pressure/transport and momentum work must be coupled
consistently before this representation is adopted in full moving evolution.

## Variable-density field component and explicit pressure bridge

`compatible_density` reconstructs actual conservative number-density means and
provides its physical quadrature adjoint. `variable_compatible_current` applies
the positive variable-density inertia by tensor products and exposes separate
local current, density and full magnetic-energy derivatives. The existing
original-equation solver consumes its matrix-free inertia without a different
physical mass, frozen-density substitute or changed tolerance. The caller still
owns continuity; this component does not integrate density or momentum.

`compatible_charge` derives weak graph current flux from the full magnetic
Hessian, preserving divergence and matching Ampere/gradient work. The deliberate
`compatible_pressure` bridge pairs the older lumped species thermodynamics with
that current. Its thermal error remains explicitly low order; nodal thermal
density and reconstructed conservative density are not silently interchanged.

Independent continuum fields, physical norms using root-product primitives,
source/array archives, original-equation replay and plots live in separate
`examples/variable_current_*` modules. VARIABLE_CURRENT_NUMERICS/VERIFICATION
separate strong field/energy convergence from the remaining entropy-rate error.
Full moving thermal compression, velocity/vorticity, stress and complete energy
flux must be derived and verified before changing the production state partition.
