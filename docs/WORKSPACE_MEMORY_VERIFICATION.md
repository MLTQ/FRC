# Workspace changes and a fully replayed axial refinement

2026-10-12 continuation from `5a60b36`. All workers are terminal and collected.
No unit tests were added or run. The complete original FRC sandbox goal and
all 12 acceptance beads remain open. These are instantaneous full coupled
signed rate calculations and independent physical audits in the original
filled periodic cylinder with impermeable radial wall and prescribed bias.
They are not time trajectories, moving-vacuum calculations or formation runs.

## Preserve arithmetic and release completed work

`ReleasedCurrentTransport` and `ReleasedRegularCurrent` add only local `del`
statements after the last use of quadrature gradient, stress, velocity, field
and current temporaries. Removing those statements leaves each method's AST
identical to the original. This preserves arithmetic expressions and call order;
the source check is supporting evidence rather than a physical verification.

`DetachedMetricWallPlasma` obtains the original fully rebound operator copies,
then removes completed-model references from cached geometry templates. Active
copies retain the actual density model throughout evaluation. The original
geometry maps, species projection and numeric approximate inverse remain.
After the complete parent derivative returns, collection releases unreachable
model/inertia cycles formed by `LinearOperator` bound methods. No state, input,
equation, physical mass/density, field basis, quadrature, boundary, residual or
budget changes. No diffusion, filtering, atmosphere or energy repair is added.

Both variants run the complete F4/64, B16/64 coarse signed fixture and the
F64/256, B256/256 fine signed fixture. Each original/new comparison requires
every saved input, coefficient, physical/face array, solver diagnostic and
q16/20 physical norm to agree bitwise. All four comparisons pass, as do their
fresh diagnostic replays. Each complete candidate calculation also has a fresh
all-array simulation replay. Different variants retain different source records;
physical equivalence is distinguished from same-source replay.

Evidence is in `released_current_*`, `detached_current_*`, and
`outputs/workspace_memory_diagnostics/workflow_equivalence.json`. Default
geometry, old calculation sources and time integrators remain unchanged.

## Measured resource limits

Own-process cumulative peak RSS includes the entire guarded physical workflow,
its independent audits, and the verifier when selected. Values vary by process;
they are not exclusive per-function allocations or guarantees for another grid.

| Complete workflow | Primary peak (GiB) | Fresh replay peak (GiB) |
|---|---:|---:|
| Previous tensor baseline, B256/256 | 7.420 | 7.150 |
| Release expired workspace, B256/256 | 6.873 | 7.430 |
| Also detach completed cache models, B256/256 | 6.446 | 7.056 |
| Detached variant, B256/288 | 7.268 | 7.408 |

The first lifetime-only result by itself would overstate dependable savings.
The cache candidate provides enough observed headroom for a modest axial
increment, not evidence that doubled or radial grids fit. Every new guard stays
untriggered, and final measured values also obey the unchanged 8 GiB own-worker
peak and 256 MiB free-disk budgets. The original larger-grid resource-interrupted
attempt remains unchanged. No process outside these workers was terminated.

Two new manifests record sharing byte-identical completed replay NPZ files after
checking committed verification records and both full hashes. Nominal duplicate
bytes total 496,702,418. All scientific bytes, paths, source archives, failed
attempts and the user's poster survive. Shared paths are immutable: never write
to either in place. Observed host free-space changes are not attributed to a
particular process or treated as a memory guarantee.

## Independent axial refinement at fixed physical reference

Both grids keep fluid F64/256, radial field 256, independent radial16,N128 mode,
epsilon 1e-4, physical masses and every original residual definition/bound.
Only axial field resolution increases, from 256 to 288. This 12.5% increment
is explicitly smaller than the pending doubled-grid study. Both signed full
solves, q16/20 volume/face audits and a fresh complete simulation replay pass.
Every saved input/rate/face array and non-runtime metric reproduces exactly.

At q20, errors below are percentages of each component's own continuum L2 norm:

| Physical rate | Axial 256 | Axial 288 |
|---|---:|---:|
| Radial magnetic field | 0.0979624% | 0.0491827% |
| Toroidal magnetic field | 0.547509% | 0.363014% |
| Axial magnetic field | 0.00994610% | 0.00518330% |
| Radial current | 5.96939% | 3.05317% |
| Azimuthal current | 0.0190466% | 0.0164420% |
| Axial current | 1.18700% | 1.11402% |
| Full current vector | 0.0878993% | 0.0817008% |

The radial-current error decreases 1.955-fold. Independent root-polynomial
differentiation of the actual magnetic rate gives radial-current errors 5.94247%
and 3.03377%. The weak/strong difference slightly increases, 0.110758% to
0.120095% of the radial reference norm. Thus the measured improvement also
appears in the magnetic response itself; it is not a corrected weak Ampere
identity. Strong curl excludes distributional face jumps. Axial weak/strong
difference remains about 0.3453%, and axial-current accuracy is still limited.

Both complete component/strong-curl reconstructions reproduce the original
physical fingerprints and fully replay. The current-error decomposition closes
within 4.23e-16 relative; largest q16/20 component-norm change is 4.76e-16.
Physical rate norms, including axis/wall regions, change by at most 1.30e-15.
The common reference L2 norms agree across the two grids within 2.23e-16 relative.
No state or spectrum was filtered. These two resolutions do not establish a
convergence order or bound the separate small-amplitude roundoff sensitivity.

The refined pair's maximum original residual is 6.63361e-11 under 1e-10.
Signed total-power defects are -1.94618e-19 W and 1.83258e-19 W. The larger
absolute defect divided by the sum of absolute fluid/current/magnetic exchange
powers is 3.111e-9. This is an instantaneous global work check; local budgets
and accumulated time-integration conservation remain independent requirements.

Full data: `outputs/detached_axial288*`, `outputs/detached_baseline_curl*`,
and `outputs/workspace_memory_diagnostics/axial_refinement_metrics.json`.
The [scientific PNG](../outputs/workspace_refinement_report_layout/verification.png)
was visually inspected; a PDF is alongside it. The first layout is preserved.

## Evidence and next action

Eleven complete verification groups comprise five signed physical calculation
replays, four original/new comparison replays, and two independent physical
curl-audit replays. They are not eleven formation experiments. Seven raw-file
comparisons cover every array in those replay records. Nine new modules have
companions; 42 new source archives, 355 preceding archives and 22 older source
closures verify. The final 156-file source union is
`7e89bbaac2a377a3abe1d2a6fab0b0641381dbd777ae04a903a5c86e8c32e918`.
Full inventory: `outputs/workspace_memory_diagnostics/checkpoint.json`.

Next reduce actual gradient/stress quadrature workspaces further if needed for
axial320+ and independent radial-field refinement under recorded resource limits.
Require full same-grid equivalence/replays for any new implementation, then
component/strong-curl, amplitude, residual and local-budget controls on each
new grid. The measured axial improvement does not settle the fixed radial/fluid
errors. Do not infer adequate accuracy from the much smaller full-current norm.

Full signed trajectories, local energy/charge/enthalpy/Poynting budgets, heated
moving exact vacuum, transport/ionization/model validity, real prescribed coils
and physical ends, converged formation studies and the experimental UI remain
open. Lab dimensions, gas, neutral fill and bias remain provisional. No3-D
stability or confinement claim; research and issue state remain local, with no
remote or publication. The original full goal is unchanged.
