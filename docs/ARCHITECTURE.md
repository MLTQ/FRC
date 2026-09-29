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
