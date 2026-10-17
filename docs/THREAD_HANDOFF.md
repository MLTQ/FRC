# Current goal continuation, 2026-10-17

The user explicitly set completion of the lab-facing virtual plasma bench as an
active goal. Retain 1 m length / 15 cm diameter / 0.1–100 mTorr and Ar/H2/D2/N2/O2.
Read [VARIANCE_TRANSFER.md](VARIANCE_TRANSFER.md) for current equations/evidence
and [PROJECT_OVERVIEW.md](PROJECT_OVERVIEW.md) for physical scope and missing work.

Actual transport-stage refinement reduces trace Ar8 maximum temperature from
7.041 to 0.391 to 0.347 eV (32/44/46 velocity cells). Full20ns heavy-only enriched
10ns-source calculation completes in971 steps with relative energy error2.84e-11.
Current remains unresolved: coarse32/64/128 differ from the existing1024-cell
numerical endpoint by66.3%/57.3%/34.2%; the reference itself is unconverged.

The first enriched5ns-source attempt22186 is EXIT1 at10ns coupled acceptance:
six retained cold Ar2 cells produce0/0 in the old sound formula after m*N
underflows. Independent100digit actual-source diagnosis is saved. A separate
full-duration run evaluates equivalent factored ratios; it does not repair heat.

Original variance5ns73344 is now EXIT0 at full20ns:972steps,3818.456s,
2415015density cells,Erel2.76935e-11. Its complete10/5ns E/J differences are
4.53748e-6/4.33645e-4, but trace Ar8 temperature differs0.161226eV RMS and
2.231827eV maximum. Final Ar8 range0.190592–5.152868eV. This remains an
accuracy limitation despite good conservation; compare the enriched cases next.

Live handles: **16195** factored-ratio heavy-enriched5ns (15ns accepted, final
interval last200steps/4.425498ns), **61204**128-bit-scaled heavy-enriched5ns
(same accepted time/latest interval clock), and **91194** corrected uniform-
spaced128-electron/32-heavy20ns case (first10ns accepts2267steps,8469cells,
Erel−4.83442e-12; second interval active). **31198 is now EXIT0** at full20ns.
All use frozen roots; exact paths are in VARIANCE_TRANSFER. Poll actual handles;
do not restart quiet jobs. Main has no remaining live source after the completed
heated neutral-collision case below. A quiet final source return can be expensive.

All six variance analytic cases57807 and enriched-subcycled analytic pair47880
are now EXIT0 at their full original durations. Independent100digit transfer
integrals exposed~0.023eV mixing heat lost by underflow in six retained cold
source cells. Binary coordinate scaling restores it without changing physical
amounts; its full20ns10ns-source case5892 EXIT0 matches original E/J to
9.43e-15/1.66e-8. The exercised arithmetic and optional coordinate shift are
now in main, after77325 completed. That full20ns e128/h32 case agrees with
the previous128-cells-per-species E/J within0.036%/0.056%, at about25% lower elapsed
time. It still differs34.2% in J from the unconverged1024 numerical reference.

Moving charged collisions now complete20ns (32422 EXIT0),971steps,89.884s,
Erel2.84193e-11. Matched reaction-only E/J differences1.56447e-6/8.96056e-4;
see ATOMIC_COLLISIONS. New ground-neutral/electron elastic force and heat use
the already authorized private BSR-500 momentum-transfer table, including recoil.
All six independent full source relaxations7161 EXIT0 show approximately
fourth-order time refinement; finest physical electron/Ar drift error1.198e-7
of initial, Te error1.708e-6eV,energy drift6.422e-8. Figure inspected. See
ELASTIC_NEUTRALS for equations, raw result paths, input/source SHA and limits.

The combined moving source also completes20ns (79419 EXIT0):971steps,103.579s,
21349cells,Erel2.84185e-11, source6539c5fcd1de3d50fc76e1ad923341d0a67808b1672f2858c5b2845d4628b8e1.
Matched charged-only E/J differences1.02593e-8/6.76392e-5; comparison99437 EXIT0
also retains all prior rows, including the full variance5ns endpoint. Neutral
collisions are explicit `--neutral-elastic-file` plus `--coulomb-log`; reaction-only
remains default. Private tables/rate curves are not bundled; no physical channel
is inferred for other gases or excited neutral targets. Bead frc-transport.7
tracks coupled accuracy and missing neutral channels; keep it open.

The original heated3.5us source now completes700/1400steps with ground-neutral
electron scattering (59917 EXIT0,84.627s), saved in
`outputs/atomic_collisions_neutral_elastic`, source
`799e2be212d00b2e4fcec7b2a767634a14b3aa21e7786d204f46fc9e444efbc8`.
Fine final Te2.502149135eV, mean Ti0.263929084eV and Tneutral0.104227706eV;
current[0,0.000605107,−29.423141409]A/m². Azimuthal current damps faster, while
the neutral/ion flow difference sustains axial current. Time-halving maxvector
difference0.188295A/m², final axial difference0.002612A/m²; maxenergy residual
3.13305e-12eV/nucleus. Replotted six-panel PNG inspected; original simulation
source and later figure provenance are separate. Readout37371 EXIT0 saves exact
flows and domain coverage without changing any history. Excited neutrals lacking
an elastic channel reach2.315% of neutrals; no ground-rate substitution.

New child**frc-transport.7.1** tracks low-energy Ar+–Ar exchange. Phelps1991
Table7 Qm covers0.1–10000eV laboratory ion energy and already includes symmetric
charge exchange. Visually inspected primary pages567–568 are at
https://srd.nist.gov/jpcrdreprint/1.555889.pdf. In the saved heated trajectory,
up to11.553% of ground ion/neutral relative-pair probability lies below0.1eV,
so a complete thermal model needs sourced low-energy behavior. No partial table,
Langevin extrapolation or double-counted charge-exchange drag is installed.
The2002 author-hosted transport paper in ELASTIC_NEUTRALS establishes the
momentum-transfer convention. Phelps2007 author notes are unfinished and do
not supply a ready validated low-energy fit. A further primary candidate is
Chicheportiche et al.2014, DOI10.1103/PhysRevE.89.063102, on ground/metastable
Ar+ transport and quantum scattering; no numerical dataset has been acquired.

The first uniform-density attempt**9098 EXIT0** completes20ns,4027steps,
110129density cells,1021.890s,Erel−4.62663e-12. Its saved initial mesh has only
102velocity cells on128uniform density cells because the constructor used
cosine targets. Comparison13317 EXIT0 gives E/J differences5.25086%/48.44737%
and electron-front difference11.25444mm from1024, versus2.44748%/34.20280%/
3.39308mm for128cosine electrons. This rejects that mixed mesh, not the intended
uniform128 velocity space. Original history and source remain intact with an
explanatory README. Initial requested counts must not substitute for saved indices.

The corrected driver explicitly supplies all initial uniform electron faces to
the mechanical constructor and records actual counts. New**91194** repeats the
full20ns from frozen `/private/tmp/frc-full-uniform-electron-20261017`, source
`1e2873c2ff590fbadb4d15c9ded1190e77579e63ee376587d78e68d0ab8e575a`, output
`reactive_material_argon_e128_full_uniform_h32_10ns`. Its initial NPZ confirms
[32,32,128] actual velocity counts and128electron spacings156.25um. Defaults
remain cosine. Fine comparison now reads actual initial counts and front shifts
for every case; the revised figure labels the original attempt102mixed.

The256-electron/32-heavy run**31198 EXIT0** completes20ns,8647steps
(4869/3778),166537density cells,3028.127s,Erel3.65072e-13. Existing fine
comparison**8942 EXIT0** reduces E/J difference from1024 to1.02547%/19.07075%
and electron-front difference to1.85810mm. Prior e128/h32 values are
2.44748%/34.20280%/3.39308mm. The32-heavy ion-front difference remains34.3293um.
The1024reference is still unconverged (512/1024 J difference8.078%); no isolated
spatial-order or laboratory accuracy claim. Updated PNG and rendered PDF are
inspected, including the corrected102mixed label and observed initial counts.

Next collect the three live refinement trajectories and compare local thermal/source-
time and initial electron-current resolution using existing comparisons. No2.5ns
launch yet. The new uniform-spacing option moves initial electron resolution
toward the measured interior error while retaining the same physical state,
count and phase cap; do not assume it improves full accuracy. Existing fine
comparison now selects the two new electron cases only after successful20ns.
Independently, source the missing low-energy ion/neutral channel and evaluate
its full physical relaxation before coupling it to the heated/moving cases.
Do not invent neutral-ion/charge-exchange data. The analytic
enriched non-subcycled attempt remains a saved failure. Defaults remain unchanged.
No unit tests, remote or publication. All parent goals/beads stay open.
Preserve the user's modified poster SHA
`a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca`.

## Previous completed continuation, 2026-10-16/17

The user wants a credible virtual plasma bench to inform construction of a real
lab. Keep the provisional 1 m length, 15 cm diameter and 0.1–100 mTorr envelope.
Available gases include H2, D2, N2, O2 and argon. Hardware remains prescribed
waveforms. The full `frc-goal` is open: component results are not validated lab
settings, FRC formation, three-dimensional stability or particle confinement.

Read [REACTION_TRANSFER.md](REACTION_TRANSFER.md) and
[MOLECULAR_GASES.md](MOLECULAR_GASES.md) for equations, evidence and limits.

## Completed this continuation

- Explained the extreme trace velocity with independent 100-digit assembly and
  an analytic counterexample. The explicit constrained-volume transfer avoids
  that specific defect, closing 5.5.4. It preserves block momentum and each fine
  cell's total energy, but local cooling/thermal accuracy remains unresolved.
  Galerkin/local_total stays default. The bounded-moment alternative is rejected
  on complete controls with 51–62% current differences.
- Affine support extension retains old independent knots and every physical
  endpoint/density cell. All six full neutral/charged analytic cases complete,
  closing the specific edge-stiffness bug 5.5.7. The MPFR family also completes:
  handle 95949 EXIT0, 2654.85 s, finest neutral/charged front errors
  1.42e-15 / 6.94e-17 m and relative energy errors 5.68e-14 / 3.27e-15.
  These analytic fixtures disclose 8-ULP cleanup; reacting argon uses zero.
- Exact polynomial force integration, optional 192-bit MPFR geometry and exact
  bracketed coordinate searches reduce observed cost without changing the
  physical equations. The integrated full 20 ns search run is EXIT0, 971 steps,
  71.96 s, with the entire saved endpoint identical to the original integrated
  MPFR trajectory. Source SHA:
  `9de32e70da00456d2d8d23fc02c4bc2cf22ceeb95afe4faf50bb6a04af03d18d`.
  The actual 227839-cell RK4 step costs 7.83 versus 12.21 s with identical output.
  The earlier two-part-key sorting experiment was slower and is not adopted.
- Corrected 32/64/128-cell reacting argon runs all complete the full 20 ns at
  10 ns source intervals. The finest run (82751 EXIT0) uses 4224 actual steps,
  84343 density cells, 1069.60 s and relative energy error 1.38e-12.
  Successive E differences decrease 6.65% to 2.70%, but J differences increase
  32.21% to 38.94%. Electron-front differences are 3.47 and 2.18 mm.
  This is not converged prediction; the plasma-phase cap also scales with cells.
- Complete matched direct-transport controls reproduce those E/J discrepancies
  with reactions/remaps disabled, dominated by electron current. Evidence:
  `outputs/reaction_direct_comparison`; handles 67617/62563 EXIT0.
  This locates the spatial limitation in existing transport under 5.5.1.
  Reuse older fine-grid/half-step evidence; do not automatically double all runs.
- H2/N2/O2 BEB ionization matches all 293 public NIST table values within printed
  rounding and completes independent quadrature refinement (scoped bead 6.1
  closed). Attributed input tables live in `data/molecular/nist_beb`.
  D2, molecular products/energy partition and complete collision channels remain
  open. No molecular gas is yet enabled as a coupled discharge model. The older
  H2 inventory switch assumes complete dissociation/ionization into protons.

## All calculations are terminal

The 5 ns-source MPFR run 76215 is terminal EXIT0: full 20 ns in 972 actual
transport steps (325/203/223/221), 5342.87 s, 2417135 density cells, minimum
width 2.68e-27 m and relative energy error 2.77e-11. The exact source is
`bfa29713ca0a32407971b23f277453a9255667e60aea688f51d39242f0ad40c1`.
Its complete archived state and endpoint readout are in
`outputs/reactive_material_argon_newborn_mpfr_5ns`. Ar9+ uses 22 independent
velocity cells, minimum physical element 0.269 mm. No old knot, density cut,
physical endpoint or represented population was removed.

The complete comparison 89248 is terminal EXIT0 and its figure was inspected.
Against 10 ns source intervals, E/J differences are 9.209e-6 / 3.792e-4 relative.
Neutral front difference is 1.731 um; electron front difference is 3.526e-11 m.
Their common-volume temperature RMS differences are 0.0008824 / 0.0001634 eV.
Ar9+ velocity/T RMS differences are 154.01 m/s / 0.01080 eV. This completes the
narrow newborn-mesh bug criteria, closing bead 5.5.8.

Thermal accuracy is NOT resolved. Ar8+ temperature differs by 0.22056 eV RMS
and 3.34554 eV maximum. Its final range is 1.10e-13 to 3.54570 eV; several other
trace ions approach 1e-13 eV. Largest last-return temperature changes across
accepted intervals are 0.00490 / 0.04966 / 7.18286 / 1.47703 eV. Total numerical
thermal change is -1.612e-11 J/m2 with zero cell-total-energy redistribution.
Good energy conservation and small net E/J sensitivity do not validate this heat.
Electron maximum pointwise T difference is 0.94810 eV despite small volume RMS;
all narrow regions remain included in both norms.

The redundant Decimal counterpart 76748 is terminal EXIT1/KeyboardInterrupt,
retired at 15 ns coupled acceptance after matching saved physical arrays and
reported fourth-interval clocks. Its termination record is explicit; it is not
counted as a full run or physical failure. All analytic, spatial, direct-control,
search-equivalence and molecular calculations are also terminal. No source root
needs to stay frozen for a live process. Original private roots remain available
and every result has a source archive. The N2 configured-entry check is an
expected unsupported-model initialization rejection, not a physical trajectory.

## Concrete next action

Fix local thermal transfer under 5.5.5 using the actual completed 5 ns trajectory,
including the large Ar8+ discrepancy. The accepted 10/15/20 ns states and frozen
source make the affected third/fourth intervals reproducible. Keep the requested
full duration, physical masses, populated species and fine density cuts when
comparing a replacement. The existing 10/5 ns postprocessor already measures
local E/J/front/velocity/T; do not build a generic verification framework.

Transport issue 5.5.1 also requires
better interior electron velocity/pressure accuracy; the complete no-reaction
controls locate that limitation independently of projection. Reuse existing
finer-grid, half-step and analytic evidence rather than automatically doubling
all runs. Full 2.5 ns-source refinement remains pending under 5.5.6; the 5 ns
case already grows to 2.42 million density cells and exposes an accuracy defect,
so a solver correction is the next milestone. Molecular physics remains bead 6.
All parent goals, including frc-goal, stay open.

The 20 ns planar n=1e14/m³ case is a numerical component and does not represent
the user's neutral-fill pressure, full chamber or demonstrated preparation.
A concise whole-project map is [PROJECT_OVERVIEW.md](PROJECT_OVERVIEW.md).

Main optional dependency: gmpy2 2.3.2, selected before import with
`FRC_GEOMETRY_BACKEND=mpfr`; Decimal remains default. Read module companions
before edits. Frozen source archives and actual steps accompany every run.
Old interrupted 5 ns attempts have explicit termination records and are not
counted as full runs. No unit tests, remote publication or viewer/poster edits.
Latest free disk is about 4.9 GiB; retain the 256 MiB reserve. Private LXCat data
stay outside Git. User poster SHA256:
`a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca`.

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
