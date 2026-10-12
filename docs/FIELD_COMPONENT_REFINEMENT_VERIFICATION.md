# Field-refinement resources and separable mass-inverse verification

2026-10-12 continuation from `09b9363`. The complete experimental FRC goal and
all 12 acceptance beads remain open. No unit tests were added or run. These
are filled, periodic, bounded-wall verification fixtures with prescribed bias;
they are not formation experiments, physical end/coil models, or demonstrations
of stability or confinement. See FIELD_COMPONENT_REFINEMENT_DESIGN.md.

## Existing fluid refinement, resolved by current component

The complete archived F32/256, B256/256 signed response now has an independent
physical current-component and strong-curl audit, with complete bitwise replay.
Both fluid resolutions use the fixed independent radial16,N128 reference,
epsilon 1e-4 and tolerance 1e-10. Relative L2 norms use each component's own
continuum norm; percentages below are physical volume errors at quadrature 20.

| Current-rate component | Fluid 32 | Fluid 64 | Error reduction |
|---|---:|---:|---:|
| Radial | 11.66354% | 5.96939% | 1.95x |
| Azimuthal | 0.135383% | 0.0190466% | 7.11x |
| Axial | 36.12045% | 1.18700% | 30.43x |

The F32 independently differentiated magnetic rate gives radial-current error
11.64947%, close to the actual current-rate error. Its axial strong-curl error
is 31.8490%, with weak/strong difference 5.36527% of the axial reference norm.
Inside-cell strong curl excludes distributional face jumps. Largest quadrature
16/20 absolute norm change is 9.72e-17. Original physical fingerprints reproduce
exactly. Evidence: `outputs/current_rate_fluid32*` and the scalar comparison in
`outputs/field_component_diagnostics/fluid_component_comparison.json`.

## Original axial refinement stopped by resource limits

Source-derived retained sparse-plus-uniform-metric estimates are 1.547 GiB for
B256/256, 3.094 GiB for B256/512 and 5.098 GiB for B512/256. These exclude full
two-dimensional LU fill, nonlinear quadrature arrays and allocator workspaces;
they are not peak-memory bounds. The small probes and exact storage formulas
are archived in `outputs/field_refinement_storage`.

The unchanged full F64/256, B256/512 initialization reached an own-process
peak RSS of 12,057,903,104 bytes (11.230 GiB) before its disk guard observed
267,878,400 free bytes. No signed input array or physical rate was produced.
Native work did not promptly yield to SIGINT. The first approval review could
not initialize because disk space was exhausted; after removing only two
agent-created temporary notes, retry approved termination. Worker 36206 exited
143 and was collected. No unrelated process or user artifact was removed.

Partial summaries intentionally retain their historical incomplete status.
The terminal interruption record is
`outputs/field_component_diagnostics/axial_interruption.json`. This is a
resource interruption, not numerical rejection, refinement success, or proof
of the system-wide cause of disk consumption. The slow original function was
not profiled; this run does not establish that global LU alone caused the peak.

## Exact separability, different floating-point evaluation

The opt-in `TensorMassGeometry` retains the original assembled electric masses,
bases, quadrature, Hessian, incidence curl and wall definitions. For each
component, M = Gr kron Gz and the row-major unknown X satisfies Gr X Gz^T = F.
Separate radial and axial LU solves replace a global two-dimensional mass LU.
All original equations and acceptance bounds still use the assembled M.
Finite-precision inverse results need not agree bitwise between algorithms.
No modes, physical masses, density, diffusion or residual budgets changed.

Complete signed F4/64, B16/64 original and tensor pilots use the independent
radial1,N64 reference. Read-only observers during fresh full pilots verify
identical fingerprints for every assembled base and strong-wall operator,
including sparse indices, coordinates, weights, masses, Hessians, curls and
boundary projections. Both observed full calculations reproduce their own
archived pilots. A fresh repeat reproduces both operator records and all
scientific arrays/metrics. This measured operator identity is coarse-grid
evidence; it is not presented as a captured fine-grid operator comparison.

The tensor fine case uses F64/256, B256/256 and the unchanged radial16,N128
reference. Both signed full solves and independent q16/20 volume/face audits
complete. A separate-process replay reproduces every saved array, physical
fingerprint and non-runtime metric. Maximum original residual is
6.49670e-11 under 1e-10. Signed total-power defects are -4.91e-20 W and
5.47e-19 W; the largest defect divided by the sum of absolute exchanged
fluid/current/magnetic powers is 8.75e-9. This is an instantaneous work check,
not a time-integrated or local conservation proof.

An independent full component reconstruction reproduces the original and new
physical fingerprints and itself replays completely. Actual old/new physical
differences at q20, normalized by the same component's continuum norm, are:

| Component | Original error | Tensor error | Actual difference |
|---|---:|---:|---:|
| Radial current rate | 0.0596938936 | 0.0596938819 | 1.04097e-5 |
| Azimuthal current rate | 0.000190465838 | 0.000190465838 | 7.13e-14 |
| Axial current rate | 0.0118700082 | 0.0118700098 | 1.36222e-6 |
| Full current-rate vector | 0.000878993012 | 0.000878993118 | 9.96358e-8 |
| Toroidal magnetic rate | 0.00547509477 | 0.00547509404 | 2.33943e-7 |
| Full electric vector | 4.56807e-6 | 4.56807e-6 | 8.25600e-14 |

Largest q16/20 absolute norm change is 3.44e-15. Fluid density and axial-velocity
differences are zero. The roughly 6% radial-current spatial error persists;
this is a verified alternative inverse evaluation, not an accuracy repair.
The original geometry remains the default and no time integrator was changed.
Full results: `outputs/tensor_mass_fine_baseline*`,
`outputs/tensor_mass_fine_comparison*`, `outputs/tensor_operator_workflow*`,
and `outputs/field_component_diagnostics/tensor_comparison_metrics.json`.

## Resource evidence and remaining work

The fine baseline and ordinary replay peaked at 7.420 and 7.150 GiB, respectively,
with neither 8 GiB / 256 MiB guard triggered. Their complete physical workflow
times were 256.34 and 272.46 seconds; these are observations, not a controlled
speedup benchmark against the earlier implementation.

The complete stage-observation replay is terminal, scientifically bitwise
verified, and recorded separately in `outputs/tensor_resource_observation`.
Its guard did not trigger. Peak RSS is 7,980,302,336 bytes (7.432 GiB). During
the first signed rate, metric construction raises cumulative peak 3.910 to
4.829 GiB, current transport raises it to 6.558 GiB, and species pressure to
6.592 GiB. The second signed electromagnetic evaluation raises the peak from
7.061 to 7.432 GiB. Thus transport is a material interval, but it is not proven
to account for the whole limiting peak. These nested high-water intervals are
not live/exclusive allocations and must not be added. Runtime/resource values
are not numerical replay metrics; every scientific result matches the baseline.

Two storage manifests record atomic hard-link sharing of byte-identical,
previously verified completed replay NPZ files. They preserve every scientific
byte and pathname, with original and repeated hashes rechecked. Nominal duplicate
file bytes total 391,609,411; observed free space varies independently. These
artifacts are immutable and must never be edited in place through either link.
No differing data, interrupted attempt, archived source or user poster was deleted.

All workers are terminal and collected. Seven new modules have companions;
37 new source archives, 318 preceding archives and 22 older source closures
are verified. The full 147-file source union is
`8ae692975fec4c2fac1345f723b13bb428c1bb1affa2c72557f527d476c8b332`.
See `outputs/field_component_diagnostics/checkpoint.json` for replay inventory,
source preservation, worker terminal states and interruption/resource records.

Next reduce measured live quadrature/workspace costs without changing parent
equations or acceptance checks, requiring full same-grid physical comparison
and complete replay for each candidate. Only then retry independent axial and
radial field refinement under recorded execution budgets. Preserve individual
current norms, actual magnetic derivatives, all strict failures and local work
defects. Full trajectories, local conservation, heated moving exact vacuum,
transport/ionization, real prescribed coils/ends, converged formation and the
experiment UI remain open. Lab dimensions, gas, fill pressure and bias remain
provisional exploration inputs. Research stays local; no remote or publication.
