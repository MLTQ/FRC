# Combined spatial refinement and complete evolution

## Verification policy

The user requests end-to-end verification only. No new unit tests are created
and no existing unit-test suite is run. Historical unit-test counts remain
historical. New evidence consists of complete signed simulation workflows,
independent physical field/conservation audits and fresh complete replays.

## Combined high-wavenumber spatial workflow

Use the frozen split_inverse_radial_rates entry point on F32/256,B256/256 and
F64/256,B256/256, where F is fluid radial elements/axial points and B is field
radial cells/axial points. Hold the independent radial16,N128 continuum
reference fixed. At each grid, execute the entire initialization, two signed
nonlinear full RHS evaluations, physical quadrature16/20 audit, initial
representation and one-sided face/current diagnostics, and frozen recording.
Repeat the complete workflow and compare every array and non-runtime metric.

The first full-array attempt was interrupted during physical reconstruction
after the host reported no space for a shell temporary file. Its traceback,
source archive and incomplete status are preserved; it produced no accepted
spatial audit. The traceback located postprocessing, not the previously
suspected matrix initialization. A concurrent trajectory replay exited after
two saved steps without a recorded cause. Its partial archive remains intact;
fresh full reruns use new destinations. No unrelated user files are deleted.

New streamed_radial_rates retains the same actual solves and saves their
inputs/results before auditing. Independent volume reconstruction uses two
radial cells per batch, retaining all axial points, quadrature16/20, physical
fields and exact face traces. Norms use row sums and math.fsum, so rounding
differences from the original audit must be measured. A complete same-mesh
baseline must reproduce every input, rate, solver diagnostic and face array
bitwise, with all norm/fingerprint differences exposed. New full workflows
and their postprocessors also receive fresh complete replays. Batching is
diagnostic only; no numerical mode, equation or accepted state changes.

The first streamed baseline completed. Its original-versus-streamed physical
fingerprints and all raw arrays/solver diagnostics agree bitwise; the maximum
absolute norm change is2.22e-16. Its rerun was intentionally stopped after
discovering repeated construction of the same axial curl basis. Both the
completed baseline and partial rerun remain frozen. cached_streamed_rates
reuses that identical basis once per interface audit, preserving all matrix
expressions. A full baseline must agree with the uncached streamed run in
every metric and fingerprint before the larger study proceeds; the cached
driver then receives its own complete fresh replay.

Both cached baseline runs and their exact replay completed. The first large
combined case also completed, taking707.845s. Before the second case returned
a result, the batch was stopped to remove another known setup cost:
WallInductionCore factors the magnetic Hessian, then UniformWallInduction
immediately replaces the result with UniformWallMetric. DirectMetricWallPlasma
constructs the identical coordinates/lift/dual-lift and used metric while
omitting only that discarded factorization. The complete coarse signed
original/new workflows, all physical diagnostics and new full replay agree
bitwise. Full large-grid equivalence is required as well.

The first large case has a separate unchanged completed-case archive, with
its parent batch explicitly marked incomplete. Exact selection replay passed;
this does not claim the interrupted second case or whole parent batch finished.
Both large grids receive fresh direct-setup runs and their own full replays.
All original sources, completed evidence and partial records remain preserved.

Both new large cases completed. The finer F64/256,B256/256 result is now
subject to smaller-amplitude and tighter-residual controls after the full
large-grid replays. These use geometry_cached_rates: within each run it reuses
only the original state-independent CompatibleGeometry for identical field
meshes/settings. Each actual fluid, state, density, wall projection, species
and induction operator is still freshly constructed. Cache scope resets for
every run/replay and records constructions/hits. Complete original/new cases
with two distinct fluid grids on one field grid agree in every array, field,
norm and solver diagnostic, and the cached batch fully replays bitwise.
This reuses unchanged factorizations; it does not replace their arithmetic.

Both full large-grid replays are now complete and bitwise verified. Finest
controls use F64/256,B256/256 at epsilon1e-5 with the original1e-10 bound,
then epsilon1e-4 with a tighter1e-11 bound. Each accepted or explicitly rejected
workflow is repeated completely; a rejected solve cannot produce a physical
rate comparison. The independent quadrature and original iterations/budgets
remain unchanged. A sequential worker interrupts itself if writable disk
headroom falls below256MiB; any such resource interruption remains incomplete,
distinct from a numerical rejection. Its guard and exact partial records are
preserved. No unrelated files are deleted and no physical tolerance relaxed.

The completed amplitude comparison also needs physical difference fields:
close error norms and raw coefficient differences answer different questions.
`current_amplitude_audit` reconstructs the saved centered J_t and E responses
using only the original sampling geometry and field mesh. It requires every
reconstructed current/E/reference/weight fingerprint to reproduce the original
full audit, then integrates their actual difference at quadrature 16 and 20.
Whole-vector and individual current-component normalizations are separate.
Per-row squared integrals and sampled maxima remain saved. Its complete
postprocessing replay checks every array and metric; it runs no plasma solver
and changes no accepted field. Fluid/magnetic differences remain separate work.

This combines the previously separate improvements to radial B and magnetic
force, then refines fluid pressure on the same field grid. The original
equations, physical masses, actual density, wall space and strict residual
bounds remain unchanged. An accepted solve is not an accepted physical error.
Preserve any rejection and its exact input; never replace a failed case with
a fitted rate. Smaller-amplitude and tighter-residual controls remain separate
requirements and must not be inferred from a mesh comparison.

## Complete time-evolution workflow

The separate split_mode_evolution wrapper runs both signs of the established
radial1,N64 independent contact benchmark through100 microseconds,60 guarded
two-stage Gauss steps on F4/32,B8/32. This is an integration control for the
new split-inverse arithmetic, not a claim of resolved radial16 dynamics.

The independent mode supplies the initial perturbation and later audit target
exp(lambda*t), never a later state, forcing or fitted update. The actual
nonlinear conservative fluid, both species, current and magnetic state evolve
together. All original stage, endpoint, positivity and induction checks remain.
Every accepted stage and endpoint is saved, hashed and replayed. A failed
trajectory retains its last accepted state and cannot pass a completion audit.

Independent quadrature16/20 reconstructs all physical fields, both pressures,
full E and J, intended growth, original total energy, particles, entropy and
positivity. The split_mode_audit wrapper evaluates E using the same actual
split-inverse RHS. Compare to the earlier same-mesh/same-duration/same-amplitude
60-step trajectory to assess the arithmetic change. The comparison evaluates
both saved states with the new RHS; it is not a bitwise cross-solver claim.
Normalize actual energy drift by the independently derived initial physical
perturbation-energy scale as well as background energy. The perturbation scale
is not asserted conserved on the nonuniform growing background.

## Physical scope and remaining gates

### Resource budget before refined evolution

The current complete trajectory inverse stores both dense fluid/magnetic
couplings, a dense magnetic operator, the fluid response and a magnetic LU
for every axial Fourier mode. At degree 5, let `nf = 36*elements - 3` and
`nm = 2*(field_nr - 1)`. Those retained complex arrays alone require
`16*nz*(3*nf*nm + 2*nm**2)` bytes. Sixty accumulated Gauss secant step/image
pairs add `2*60*2*(nf + nm)*nz*8` bytes. These are source-derived storage
estimates, not measured peak memory or numerical verification.

| Fluid / field grid | Dense modal arrays | Full secant history |
|---|---:|---:|
| F8/64, B32/64 | 0.0579 GiB | 0.0397 GiB |
| F16/128, B64/128 | 0.474 GiB | 0.160 GiB |
| F64/256, B256/256 | 15.414 GiB | 1.287 GiB |

Sparse LU, field masses, nonlinear operators, Krylov and least-squares
workspaces, temporary construction and physical audits are additional.
The finest instantaneous-rate run therefore does not establish feasibility
of a full trajectory on that grid. Before launching it, budget the complete
solver and audit footprint. A next complete intermediate-grid comparison
can use F8/64, B32/64 against the preserved F4/64, B16/64 trajectories,
keeping the radial1 reference and all original physical inputs fixed.
This candidate grid is not presumed accurate: both signs, a complete transit,
independent field/conservation audits, full replay and time controls remain
required. Preserve any failed or growing response. The exact formulas and
source hashes are in trajectory_storage_estimate.json in the joint diagnostics.

Both workflows retain the filled periodic argon verification fixture, fixed
impermeable radial wall, prescribed constant bias and physical finite masses.
See RADIAL_MODE_DESIGN.md for the full original parent equations and parameters.
No explicit numerical diffusion, resistivity, filter, atmosphere, artificial
mass or energy repair is added. Spatial approximation and unchanged Gauss time
discretization still have measured or unresolved dispersion/diffusion errors;
unit-modulus linear Gauss stability is not a phase-accuracy guarantee.

Full mesh/time/residual convergence and local current/enthalpy/Poynting face
budgets remain necessary. Heated moving exact vacuum, transport/ionization,
physical coil/end conditions and converged FRC formation remain unfinished.
The user's provisional1m quartz tube and0.1-100mTorr neutral fill are not these
benchmark parameters. Neither closed flux nor an axisymmetric calculation
establishes three-dimensional stability or particle confinement.
