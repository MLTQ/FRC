# FRC simulation — current continuation, 2026-10-16

Read [REACTION_TRANSFER.md](REACTION_TRANSFER.md) for numerical evidence and
[MOLECULAR_GASES.md](MOLECULAR_GASES.md) for the user's gas requirements. The full
`frc-goal` remains open. This is an exploratory simulation project, not yet a
validated basis for constructing a formation experiment.

The intended lab envelope remains approximately 1 m long, 15 cm in diameter,
and 0.1–100 mTorr. Available gases include H2, D2, N2, O2 and argon. Argon is the
present numerical baseline; molecular gases require actual channel and transport
data, not mass substitutions. Hardware remains prescribed waveforms.

## Completed numerical work

- Independent 100-digit assembly reproduces the extreme trace-ion velocity.
  Density-weighted approximation across steep contrast is the cause, rather
  than the linear solve. The explicit `constrained_velocity` alternative
  preserves block momentum and each fine cell's total energy. It avoids that
  extreme speed, closing the specific bug 5.5.4. Local thermal accuracy remains unresolved. The default is
  still Galerkin/local_total; the bounded alternative is rejected for 51–62%
  current error in complete zero-source controls.
- All six full neutral/charged analytic expansion cases pass with affine edge
  extension. Finest neutral front error is 1.42e-15 m, energy error 5.68e-14
  relative. A separate Galerkin extension control passes. Bead 5.5.7 is closed
  for this specific support-edge defect, not for general reaction accuracy.
- Exact polynomial electric-force integration retains every density cut and
  the same integral. A complete 20 ns comparison takes 144 versus 233 seconds;
  relative E/J differences are 1.86e-14 / 1.70e-8.
- Complete 32/64-cell argon runs at 10 ns source intervals conserve energy to
  2.81e-11 / 5.60e-12, but differ in E/J by 6.65% / 32.21%. The electron front
  differs by 3.47 mm. Spatial refinement also tightens the plasma-phase cap.
  These are not converged physical predictions.
- Combined extension/polynomial/constrained transfer completes the full 20 ns
  at 10 ns source intervals in 971 steps. A separate 192-bit MPFR geometry
  implementation completes the identical case in 75 versus 136 seconds, with
  identical saved physical arrays and only approximately 1e-50 m coordinate
  differences. Its full analytic family is still running. The backend is now
  available explicitly in main; Decimal remains the default. The integrated
  full 20 ns replay matches the isolated prototype exactly.

## Current solver and live calculations

The distinct newborn-mesh defect is bead 5.5.8: nominal neighboring targets can
select almost coincident density faces. The old Ar9+ mechanical interval was
1.27e-9 m despite nominal spacing 4.18e-5 m, then compressed to approximately
1e-12 m. An explicit target-relative spacing rule now exists in main and in
frozen source copies. It affects only newly created velocity spaces, never
existing knots or density cells. Its default is zero; the experiment uses 0.5.
Actual newborn dimension can decrease and is recorded. Local accuracy remains
an acceptance requirement.

Live handles (inspect before changing their imported source):

- 76748, `/private/tmp/frc-newborn-mesh-20261016`, output
  `reactive_material_argon_newborn_5ns`: 15 ns coupled time accepted of the full
  20 ns target, 227,839 density cells; fourth interval is active.
- 95949, `/private/tmp/frc-mpfr-20261016`, output `precise_material_mpfr`:
  six full analytic cases; neutral 16/32 cells complete, 64 cells in progress.
- 76215, the same MPFR source, output `reactive_material_argon_newborn_mpfr_5ns`:
  15 ns of 20 ns accepted; saved physical arrays match the Decimal run exactly
  through that endpoint. Both users of this source root must remain frozen.
  The original isolated dependency is `/private/tmp/frc-gmp-runtime`.
- 6777, main repository, output `reactive_material_argon_integrated_mpfr_n64_10ns`:
  corrected 64-cell full 20 ns target with 10 ns source intervals; first 10 ns
  accepted, second interval in progress. Keep main imported Python sources frozen.

Main's optional dependency gmpy2 2.3.2 is now installed in `.venv` and pinned by
`plasma-bench[precision]`. Select `FRC_GEOMETRY_BACKEND=mpfr` before process start;
192-bit geometry is explicit. Current source SHA:
`6148d66b931a92019ea5e14fa7c74a6ea16e945ef16e1dd67f2f72c97f369d3a`.
The newborn-mesh analytic control 93612 is terminal EXIT0: both complete
16-cell neutral and charged cases, with inspected plots. It is not live.

The newborn MPFR 10 ns-source calculation is terminal EXIT0: full 20 ns,
971 steps, 72.34 seconds, relative energy error 2.84e-11. Compare its local
fields/fronts/temperatures with the original selection using the existing
`reactive_material_argon_constrained_comparison/compare.py` study.

Old 5 ns jobs 77785, 19305 and 42285 were deliberately interrupted at 10 ns
accepted coupled time after diagnosis and replacement. Their `termination.json`
files record actual KeyboardInterrupt exit 1. None is a full-duration success
or a terminal physical failure. Diagnostic prefixes are explicitly distinct
from accepted coupled intervals.

## Next concrete work

Collect the full newborn 5 ns and MPFR analytic endpoints. Measure complete
10/5 ns E/J, front, velocity and temperature sensitivity with the existing
comparison. Collect the corrected 64-cell spatial comparison and the full analytic
results; preserve exact geometry and every population. Keep beads 5.5.5/.6/.8
open until their actual criteria are met. Molecular physics is `frc-transport.6`.
Return to physical collision, heat, molecular and formation coupling after the
specific solver corrections, without expanding a generic verification framework.

Read the companions for precise_material, precise_reactions,
constrained_velocity_projection, velocity_extension, velocity_selection and
subcycled_material before editing. Every output has its frozen source archive.
No unit tests, publication or viewer/poster changes. Latest free disk: 5.8 GiB;
keep the 256 MiB reserve. There is no Git remote. Runtime and preserved poster
hash are in the historical material below; current status above takes priority.

## Historical handoff before this continuation

This is the current starting point for the new session. The user requested a
bookmark and next-round beads; no new simulation was launched during wrap-up.
The overarching **frc-goal remains open**. The October 4 handoff is historical
and recoverable from Git; its instructions and old test counts are superseded.

## Start here

1. Read `AGENTS.md`, the latest section of `GOAL_STATUS.md`, and this document.
2. Run `bd --no-daemon show frc-transport.5.5.4`. It is the next concrete task.
3. Read the companions for `precise_reactions`, `precise_material`,
   `mass_scaling`, `projection_energy`, and `reactive_material` before editing.
4. Use the saved failure in
   `outputs/reactive_material_argon_retained_knots_10ns/rejected_source.npz`.
   It includes exact `velocity_knots_species*_decimal_m`; do not reconstruct
   the old cosine velocity mesh and call that a replay of the inherited mesh.

Workspace: `/Users/max/Code/FRC`; no configured Git remote. Save work and beads
locally. Research base is `ab3103eec`; the following handoff checkpoint preserves
the candidate source, terminal run and snapshot evidence described below. Use
`git log -3 --oneline` to find that checkpoint rather than treating the base as
current HEAD. All referenced output paths are relative to this repository.

## Direction and user preferences

Build a physically credible, interactive FRC experiment sandbox, preserving
[the full original objective](ORIGINAL_PROJECT_REQUEST.md). The provisional lab
envelope is a roughly 1 m quartz tube, 10–20 cm diameter, selectable gas,
0.1–100 mTorr neutral fill, and adjustable Helmholtz-like bias. These are inputs
to model selection, not validated experiment settings.

Prioritize simulation development and complete end-to-end physical comparisons.
**Do not add unit tests or run the existing unit suite.** Avoid expanding
verification/reporting frameworks. Use independent mathematical or experimental
baselines to evaluate concrete changes. No subagent delegation is requested.

Retain fine resolution. Physical-width coarsening is off; the new charged driver
also defaults to zero ULP cleanup. Arbitrary run limits do not demonstrate
hardware exhaustion or justify coarsening. Preserve exact vacuum, physical
masses, all populated species and support endpoints. No temperature/mass floors,
heat clipping, discarded traces, global energy resets, hidden timestep retries
or shortened-duration success claims. Numerical mixing must be measured.

Hardware stays prescribed waveforms. Current planar Euler–Poisson/reaction
components are not FRC formation or lab recommendations. Closed flux does not
prove three-dimensional stability or particle confinement. Physical collisions,
transport, cylindrical geometry, coils/ends and converged formation remain open.

## What is established

- Diagonal mass scaling solves the same physical equation `Mv=b`, including
  acceleration, without changing mass or resolution. Earlier failed-state
  solutions agree with independent 100-digit elimination within 2.65e-15 in
  normalized maximum velocity. The current local error still needs diagnosis;
  a well-conditioned scaled matrix alone does not guarantee an accurate load.
- Fifty-digit position/difference arithmetic preserves distinct moving supports.
  All six full neutral/charged affine cases completed in
  `outputs/precise_material_expansion`. Finest neutral front error 4.79e-16 m;
  finest charged relative energy error 3.27e-15. Geometry-specific bead
  `frc-transport.5.5.3` is closed. This does not close the coupling or FRC goal.
- Retaining existing moving velocity knots eliminates a major transfer loss.
  Full matched zero-source 20 ns controls in
  `outputs/reaction_remap_retained_knots` reduce current differences from
  55.020%/65.369% to 0.01021%/0.02575% for 10/5 ns source intervals. Electric-field
  differences fall to 3.14214e-7/3.09182e-7 relative. Residuals include transfer
  and step grouping; no formal order claim. Every density cut is retained and
  actual velocity-cell counts may grow.
- The earlier `reactive_material_argon_balanced_10ns` case completed 20 ns,
  815 transport steps, relative energy change 4.24e-11, but used rebuilt
  velocity meshes. Its 5 ns and 64-cell refinements failed. It is not converged
  adoption and does not supersede the accurate-current transfer work.

## Latest terminal run and exact defect

`outputs/reactive_material_argon_retained_knots_10ns` is terminal **EXIT 1**
(handle 17499), elapsed 230.838 s. Source SHA256:
`552213d640f98cfcf49d5679fcdc6535050e940b31a62f42dda4dc8ff66f46ef`.
It requested 20 ns with 32 initial velocity cells, 10 ns source intervals,
retained knots, precise geometry and zero cleanup. The first 10 ns is accepted
(527 transport steps). All 460 transport steps of the second 10 ns interval
completed, but its final source return rejected negative thermal energy.
**Accepted coupled time is 10 ns, not 20 ns.** The retained-knots 5 ns attempt
also failed; it predates saved-knot metadata and requires evolution replay for
an exact diagnosis of its rejected projection.

The latest failure is in **inherited neutral species 0**, not only newborn
ions. Source temperature 0.00775447 eV and velocity 2377.73–2408.39 m/s become
projected velocity 2647.43–2656.04 m/s. With `U=E-K`, 16 fine cells are negative,
minimum temperature −0.17233045 eV; balanced mass condition is 6.1192.
Local refinement at all negative-cell endpoints increases velocity cells from
38 to 52, but leaves four negative cells and Tmin −0.02318735 eV even when all
failed-cell endpoints are already knots. This was a snapshot experiment;
no adaptive production refinement has been implemented.

Read `summary.json`, `exact_projection_diagnosis.json`,
`local_refinement_diagnosis.json` and `local_refinement_knots.npz` in that output
directory. `failed_transport_stages.npz` is a split-stage record, not an accepted
coupled endpoint. Raw initial/accepted state, actual transport steps, exact
rejected source and original trajectory source archive are all retained.

## Candidate energy transfer: saved, not adopted

`energy_projection='local_total'` remains the default. The explicit research
option `'projected_frame'` uses `U_c=E_c-W_c+K_c`, where `W_c=integral(p·v)`;
local total-energy redistribution is `2K_c-W_c`. Its sum is `vᵀ(Mv-b)` and
vanishes for an exact Galerkin solve. With constant source density within a
fine cell, thermal energy equals source heat plus unresolved velocity variance.
This is numerical mixing, not physical atomic heating. General density
variation and floating-point cancellation still require scrutiny.

On the exact rejected snapshot, this option yields positive temperature and
relative total-energy difference 3.03e-15, but **Ar9+ reaches 2.882045e6 eV**
from a source near 0.199 eV. Source velocity 275.6–313.6 m/s becomes target
endpoints −51949.1 and 7941858.5 m/s, despite balanced matrix condition 3.0434.
The worst cell has relative number 1.1691119e-41 with species scale
4.4727663e-270 /m². Diagnose arithmetic/load scaling and velocity approximation;
do not discard this population or accept global conservation as local accuracy.

Evidence: `projected_frame_diagnosis.json` and
`projected_frame_local_errors.json`. Global thermal-gain diagnostic is
−5.61e-23 J/m² (tiny signed roundoff), local energy-redistribution L1 is
1.92e-22 J/m², and maximum temperature change is 2.882045e6 eV. The maximum
exposes what small global totals hide. The candidate has **no complete E2E run**,
and the example driver has no dedicated CLI switch for it yet. Do not spend
this next session merely wiring it into the UI or accepting it because it runs.

Candidate source is preserved separately in `projection_candidate_sources.tar.gz`
with `projection_candidate_provenance.json`; it overlays the trajectory's
archived source. Inline diagnosis scripts were not saved; exact input snapshots,
knots and measured outputs are saved. Ledger additions are diagnostics only,
not physical energy reservoirs or compensating corrections.

## Next-round beads

| Bead | Next action / dependency |
|---|---|
| **frc-transport.5.5.4** | Start here: explain and correct the extreme trace-ion velocity projection error using actual saved data and an independent moment/projection comparison. |
| **frc-transport.5.5.5** | After .4, resolve locally accurate, admissible momentum/energy return; the neutral failure also remains. Derive and measure numerical mixing and energy redistribution. |
| **frc-transport.5.5.6** | After .4/.5, complete original 20 ns charged trajectories with source intervals 10/5/2.5 ns and spatial refinement. Measure J/E/front/temperature and budgets, varying one approximation at a time where possible. |

These are children of `frc-transport.5.5`; .5 depends on .4 and .6 on both.
`frc-transport.5.5.2` remains the in-progress reaction/motion coupling umbrella;
`.5.5.1` retains front/current accuracy work. Keep parents and `frc-goal` open.
After these concrete solver fixes, return to physical moving collision coupling
under existing `frc-transport.3`; do not create a new verification framework.

## Execution and preserved state

No relevant numerical process was found during handoff inspection. Older LIVE
entries in append-only status/history documents are historical. The long
16000-step formation pulse 17262 stopped at 1994 accepted steps, 436.1875 ns of
3.5 µs, due to its disk-space guard; it is not a completed refinement. Earlier
fine-vector and optional coarsened runs are also terminal. Do not duplicate them
based on stale progress notes. Recheck storage before large runs: about 3.5 GiB
was available at handoff, and the existing 256 MiB reserve remains in force.

Use `.venv/bin/python` with `MPLBACKEND=Agg`,
`MPLCONFIGDIR=/private/tmp/frc-mpl-cache`, `XDG_CACHE_HOME=/private/tmp/frc-plot-cache`,
`OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1`.
Use existing `examples.reactive_material_subcycling`, `reaction_remap_control`
and their companion docs for subsequent concrete E2E runs; choose a new output
directory and preserve source provenance. No unit tests were run in wrap-up.

The detailed 3-D viewer at `http://127.0.0.1:8790/20mt/` remains based on saved
formation states. No visualization work is pending, and these planar component
results must not be presented as updated FRC formation. Preserve the user's
modified `outputs/frc_showcase/poster.png`, SHA256
`a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca`.
It is excluded from the checkpoint. Private LXCat inputs stay outside Git;
prior acceptance was for local research, not redistribution. No publication.
