# Active continuation: trajectory rejection

## Active trajectory rejection replay (2026-10-12, after c912fdc)

Latest completed work: independent density localization and its fresh full
replay reproduce all216 original density/error fingerprints and all3 derived
raw files bitwise. From13.333us onward more than99.7% of the refined squared
density error is in low axial cell modes abs(k)<=4. Final error spans radius;
first-interface jump is2.27545% in background-normalized surface L2, sampled
maximum4.69502%. Localization does not identify a cause or prove local flux
failure. Initial tiny-error partitions are quadrature-sensitive. Read the new
TRAJECTORY_REFINEMENT_VERIFICATION section and inspected scientific report in
outputs/trajectory_rejection_report. No extra physical trajectory was launched.

After the complete failed-trajectory replay, the proposed next diagnostic is
the full unchanged spatial RHS along the saved early error direction, with
amplitude and independent-mode controls. Read TRAJECTORY_TANGENT_PROBE_DESIGN.md;
this probe is not implemented or executed. Opposite-sign, same-grid arithmetic,
time and complete physical/conservation controls remain required.

The refined F8/64,B32/64 negative100us/60-step attempt rejected step18 after
17 accepted steps (28.333333us). Newton GMRES failed with info=5. The original
stage1e-12/relative1e-10 bounds remain. Session36840/PID89880 is terminal,
exit0 collected through the driver's recorded numerical-failure path. Its
4.439148GiB peak/390.625MiB minimum free disk pass the8GiB/256MiB guards.
This is a failed full-duration gate, not a resource interruption or success.

Independent density-only observations through all17 accepted steps, quadrature
16/20, fully replay every time-series byte and physical fingerprint. At the
last endpoint, the refined density departs16.0233% at a stored node; its L2
error is790.258 times the growing independent linear density signal. Older
medium/coarse controls have relative errors5.700e-5/3.898e-4 at the same time.
Integral quadrature agreement is close; sampled maxima are less resolved.
Single-sign density evidence does not establish causality or full physical
accuracy. Read TRAJECTORY_REFINEMENT_VERIFICATION.md and the companion for
examples.trajectory_density_prefix. Stored terminal initial/final state and
magnetic arrays match the initial/last accepted arrays in dtype, shape and bytes.

LIVE: parent unified session40737/PID93211, replay child PID93212. Parent record:
outputs/trajectory_rejection_campaign/summary.json. It runs only the full
negative rejection replay, with original initial data, equations and strict
checks. Latest snapshot:12 accepted replay steps, not completed verification.
The one-job plan is outputs/trajectory_refinement_diagnostics/rejection_replay_plan.json.
The old seven-job plan is preserved but was NOT launched. The controller and
new density observer are now used/frozen; do not edit their Python sources.

NEXT poll the actual parent40737, without restarting or launching another heavy
job. Collect its terminal exit and verify all initial/final/state/stage files,
non-runtime solver/failure metrics and exact failed step before claiming a
reproduced rejection. Reproduced failure still fails100us. Inspect the result
before selecting opposite-sign, same-grid arithmetic and120-step controls;
retain full requested durations. Do not launch the obsolete seven-job plan as
is, because its first replay is already active. Original preserved medium
controls use older induction arithmetic, so a same-grid split control is needed.
Full signed field/current/E, pressure, conservation/history/energy audits and
mesh/time controls remain open. Resource-budget each later job and preserve
all failures; full remaining campaign storage is not presumed available.

No unit tests. The full goal/all12 acceptance beads remain open. Completed
source and evidence form a local progress checkpoint; evolving live replay
outputs remain outside that completed evidence set.
All168 old sources,251 prior evidence files,504 prior source archives and the
user's poster were rechecked unchanged. Active173-file source:
8dc9d2e9557b53a27351c973e71762c0577837e4b0024bab1c718f7b209066d8.
See localization_preservation.json and active_source_v2 for the new closure.
See primary_terminal_preservation.json, primary_exit.json and live_progress.json
under outputs/trajectory_refinement_diagnostics. Earlier all-workers-terminal
statements below are historical. This is progress, not a blocked or complete goal.

---

# Continue from complete refined amplitude controls

All workers are terminal and collected. No session needs resuming. No unit
tests were added or run; retain the user's end-to-end-only preference.
The full FRC goal and all12 acceptance beads remain open. This continuation
from cd66ec5 made verified progress.

Read the newest GOAL_STATUS, REFINED_CONTROL_DESIGN/VERIFICATION and
outputs/refined_control_diagnostics/checkpoint.json. Scientific PNG/PDF:
outputs/refined_control_report; PNG visually inspected. The old live snapshots
under refined_control_diagnostics are historical; physical_completion.json and
the final checkpoint supersede them. Do not restart those finished jobs.

## Verified results and limits

- Four accepted amplitudes: epsilon1e-5/1e-3 on B256/320 and B320/256, fixed
  F64/256 and independent radial16/N128 reference, original tolerance1e-10.
  Every full signed workflow and fresh separate-process replay passes.
- Both grids at epsilon1e-4/tolerance1e-11 strictly reject on minus after all
  five outer attempts; final dual-Faraday6.14598e-11/8.29280e-11. Complete
  rejection replays match all input arrays and metadata. Saved rejected inputs
  equal the baseline in shape/dtype/bytes, not a post-call immutability proof.
- Four component, two spectrum and two magnetic/current-curl audits replay.
  Axial-grid Jr error1.76199% becomes4.17506% for epsilon1e-5; strong curl
 1.73644% becomes3.19274%. Radial-grid Jr5.90438% becomes6.70354%. Actual
  amplitude differences3.83376%/3.19588% are hidden by small total-J errors.
  Larger amplitudes slightly worsen radial/axial current and total-J errors;
  azimuthal current barely changes. No amplitude-convergence claim.
- Above cell mode|k|64, Jr error grows9.9952x/9.7616x as epsilon shrinks10x.
  That band holds93.7982%/87.2576% of squared amplitude difference. This fits
  numerical noise amplified by differencing but establishes no cause. Lower-
  band errors persist. Cell spectra fold aliases; strong curl omits face jumps.
- All240 spectrum/component and48 current/component-curl integrals match.
  Volume component quadrature differences<=3.440e-15; band differences reach
 4.968e-5. Parseval/curl decomposition defects<=3.150e-16. Instantaneous
  work defect<=3.983e-8/exchanged powers; local and time budgets remain open.
- Maximum physical peak7.574GiB, minimum physical disk605.469MiB satisfy all
 8GiB/256MiB guards including final extrema. Audit resource checks occur at
  stage boundaries only. All shared completed replay paths remain immutable.
- Fourteen replay groups /14 raw comparisons;58 new archives,446 preceding
  archives,164 preceding sources,228 prior evidence files and22 active older
  source closures verify. Four new modules have companions. Final168-file SHA:
 5deb5dc6c4a3229d60d13d33b2a182716aadd61fcd1fcfa1fa5264ef6dff4966.

## Concrete next action: complete trajectory refinement

1. Read JOINT_E2E_DESIGN.md, PHYSICAL_INTERFACE_VERIFICATION.md and companions
   for matched_mode_fixture, matched_mode_evolution, split_mode_evolution and
   split_mode_audit. Used Python sources are frozen; add a new guarded opt-in
   orchestration module if necessary, with explicit dynamic source seeds.
2. Use the unchanged split_mode_evolution full nonlinear protocol on
   elements8,points64,field_nr32,field_nz64, retaining the existing radial1/N64
   reference and original physical settings. Run both signed epsilon1e-4
   perturbations for100us with60 and120 Gauss steps, original stage/endpoint/
   induction bounds and complete saved stage/state records. Preserve failures;
   a partial run cannot pass the full-duration gate. Run one large job at once.
3. Compare to preserved F4/64,B16/64 trajectories
   outputs/matched_mode_guarded_ne4_z64_t60_minus and _plus, keeping physics,
   duration, reference and amplitude fixed. A solver-arithmetic difference must
   be explicitly audited, not called a bitwise cross-implementation match.
   Fresh full replays of each new trajectory, independent q16/20 field/pressure/
   current/electric/energy/particle/entropy audits and time differences are
   required. The older4.31611% final-current error and fast near-grid branch
   remain unresolved; refining is not presumed to repair them.
4. Source-derived retained dense modal storage0.0579GiB and60-step secant
   history0.0397GiB are lower estimates, not peak limits. Keep actual8GiB RSS
   and256MiB disk guards and final-extrema checks; budget all saved states,
   replays and Git evidence. The finest instantaneous grid is not a feasible
   trajectory simply because its signed rates fit memory.
5. Keep refined-grid noise/strict1e-11 rejection and joint B320/320 resource/
   accuracy controls open. Any compensated arithmetic or expanded electric
   representation requires unchanged real equations/adjoints/bounds and full
   physical controls. No filters, floors, artificial mass or repairs. Continue
   local budgets, heated moving exact vacuum, transport/ionization/validity,
   prescribed real coils/ends, converged formation and experimental UI.

Preserve every archive and the user-modified outputs/frc_showcase/poster.png,
SHA a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca.
No remote or publication; local issue export only. Provisional tube/gas/fill/
bias inputs are not validated laboratory recommendations. No3-D stability or
particle-confinement claim. Older handoffs below are historical.

---

# Continue from independently refined field directions

All workers are terminal and collected. No process or session needs resuming.
This continuation from a67d240 made verified progress. No unit tests were added
or run; continue with complete end-to-end physical workflows, independent
comparisons, convergence/conservation measurements and fresh full replays.
The full original FRC sandbox goal and all 12 acceptance beads remain open.

Read newest GOAL_STATUS, SCHEDULED_WORKSPACE_DESIGN/VERIFICATION and
outputs/scheduled_workspace_diagnostics/checkpoint.json. Scientific PNG/PDF:
outputs/directional_workspace_report; PNG visually inspected.

## Verified evidence

- Scheduling alone preserves coarse/fine physical results bitwise but peaks
  6.481/6.476 GiB. A read-only full replay identifies a second transport peak
  during full stress construction. ComponentDriftStress supplies original
  ((m_reduced/e**2)*J_i)*J_j/n entries to the unchanged stress adjoint, retaining
  arithmetic association and independent off-diagonal entries. Complete
  original/new coarse/fine comparisons and all candidate replays are bitwise.
  Fine primary/replay peak 5.342/5.516 GiB; defaults/time integrators unchanged.
- Fixed F64/256, independent radial16/N128 reference, epsilon1e-4, bound1e-10.
  New B256/320 axial and B320/256 radial signed workflows and full fresh replays
  pass every original check and independent q16/20 physical/face audit.
- Axial refinement Jr errors: 5.96939% (256), 3.05317% (288), 1.76199% (320).
  Independent strong-curl Jr at320 is1.73644%; weak/strong difference0.131541%
  slightly grows. Br error falls to0.0266691%. Modest-grid log slopes5.692/5.218
  do not prove asymptotic convergence or resolve amplitude/cancellation error.
- Radial320 instead gives Jz0.734849% versus1.18700%, strong-curl Jz0.599842%
  versus1.00320%, and weak/strong Jz0.210047% versus0.345296%. Jr stays5.90438%.
  Radial-velocity error improves0.183143% to0.0983390%. Total-J0.0552027%
  still hides weaker components. Strong curl excludes face distributions.
- New maximum residuals6.82535e-11/8.96429e-11 stay below1e-10. Largest
  instantaneous power defect4.343e-9/exchanged powers; local/time budgets open.
  New q16/20 component differences<=2.64e-16; physical rate norms<=1.39e-15.
- Axial primary/replay peak6.901/7.162 GiB; radial7.161/7.037 GiB. An earlier
  axial replay hit the disk guard (262144000 free bytes), exit130 collected.
  Partial physical/source/control records are preserved; no solver rejection.
  Fresh standalone replay succeeds. Final successful guards all obey8GiB and
  256MiB limits; do not run a large solver alongside a diagnostic audit.
- Twelve replay groups: six complete signed simulations, four implementation
  comparisons, two independent curl audits. All8 raw comparisons,49 new archives,
  397 preceding archives and22 older active source closures verify. Eight new
  modules have companions. Final164-file source SHA:
  a7022b5b4d1a8f5b598c431e71f6865d51c178ca61eed43a33ef230cfaeeabf5.

## Concrete next action: refined controls, then joint accuracy

1. Read component_stress_transport, guarded_component_stress, current_rate_curl_audit
   companions and previous amplitude/curl/residual verification. All used Python
   files are frozen. Existing guarded_component_stress accepts nr/nz but fixes
   epsilon1e-4 and tolerance1e-10. Add a NEW explicitly sourced opt-in protocol
   for epsilon1e-3/1e-5 and tighter residual controls on B256/320 and B320/256.
   Keep both signed full rates, exact actual inputs, every original residual,
   q16/20 physical/face norms, component differences and full fresh replays.
   The preceding 1e-11 rejection and raw-curl arithmetic evidence remain valid;
   do not loosen a bound or assume refinement repairs that cancellation.
2. Preserve component and actual magnetic-curl norms. Separate amplitude
   differences from continuum error; the small total-J norm can conceal Jr/Jz.
   Keep point/axis/wall and initial face diagnostics distinct from evolved
   local charge/enthalpy/Poynting balances. No filters, atmosphere or repairs.
3. Budget combined B320/320 and further independent refinement only after
   measured resource headroom. Radial field count must subdivide64 fluid cells.
   The unchanged8GiB/256MiB guard must include final observations. A combined
   run's feasibility is unproved. If storage is limiting, profile the component
   variant in a required full replay before choosing another workspace change;
   original UniformWallMetric modal factors grow roughly as nr^2*nz.
4. Continue full signed trajectory mesh/time/residual convergence and local
   budgets. The radial1 F8/64,B32/64 full100us control remains pending. Rate-only
   fidelity does not establish time stability or permit formation adoption.
   Then resume heated moving exact vacuum, transport/ionization/model validity,
   real prescribed coil and physical-end boundaries, converged formation and UI.

## Data and execution

New physical source:000fd344b524568a655485308fb72d754aa8c144d26cf3a8e04c4d947ffb53ef.
Fine baseline: outputs/component_stress_fine; axial: component_axial320;
radial: component_radial320. Successful axial replay is suffixed _replay_serial,
not the preserved interrupted _replay. Each curl audit has _curl, _curl_replay
and _curl_verification directories. Diagnostics include complete source/replay
inventories, directional metrics, interruption and verified storage manifests.

All byte-identical shared NPZ paths are immutable; never overwrite any in place.
A cross-study inventory rehashed1776 completed paths before/after sharing852
additional duplicate inodes. No scientific content/path or user file was lost.
Do not delete partial attempts or retry into an existing destination. All old
source/evidence and the user's modified outputs/frc_showcase/poster.png survive.
No remote or publication. Provisional lab dimensions, gas, neutral fill and bias
remain exploration inputs; no3-D stability or particle-confinement claim.
Older handoffs follow as historical records, not active worker instructions.

---

# Continue from fully replayed axial refinement and workspace verification

All workers are terminal and collected. No handle needs resuming. No unit tests
were added or run; continue with full end-to-end physical workflows, independent
comparisons, convergence, conservation and complete replays. The full original
FRC goal and all12 acceptance beads remain open. This turn and5a60b36 both
made verified progress.

Read newest GOAL_STATUS, WORKSPACE_MEMORY_DESIGN/VERIFICATION, and
outputs/workspace_memory_diagnostics/checkpoint.json. Scientific PNG/PDF:
outputs/workspace_refinement_report_layout; PNG visually inspected.

## Verified evidence

- Releasing expired current-workspace temporaries and detaching completed
  density models from geometry cache templates preserves every original coarse/
  fine input, rate, face array, solver metric and independent q16/20 norm bitwise.
  Four original/new whole-workflow comparisons and fresh candidate replays pass.
  Original functions, default geometry and time integrators remain unchanged.
- Finest baseline remains F64/256,B256/256 on independent radial16,N128,
  epsilon1e-4,tolerance1e-10. The detached-cache baseline and full replay peak
  6.446/7.056GiB; lifetime-only6.873/7.430GiB shows resource variability.
- New F64/256,B256/288 full signed case and exact complete replay pass all
  original checks, worst6.63361e-11 under1e-10. Peaks7.268/7.408GiB obey the
  unchanged8GiB/256MiB guard; no guard fires. Do not infer doubled-grid feasibility.
- Independent current-component and actual magnetic-rate curl audits reproduce
  original fingerprints and replay fully. Jr error5.96939%→3.05317%; strong-curl
  Jr5.94247%→3.03377%. Weak/strong difference0.110758%→0.120095%. Thus the
  actual magnetic response improves, without a repaired weak Ampere identity.
  Jz remains1.11402%; full-J error0.0817008% hides weak components. Br error
  halves to0.0491827%; Btheta improves to0.363014%. Two grids are not an order.
- Refined instantaneous power defect<=3.111e-9 of exchanged powers. Largest
  q16/20 current-norm change4.76e-16 absolute; physical rate norms1.30e-15.
  Strong curl excludes distributional jumps; local/time budgets remain open.
- All11 replay groups and7 raw-file comparisons pass;9 new modules have
  companions;42 new source archives,355 earlier archives and22 older closures
  verify. Final156-file source SHA7e89bbaac2a377a3abe1d2a6fab0b0641381dbd777ae04a903a5c86e8c32e918.

## Concrete next action: memory and directional accuracy

1. Inspect released_current_transport.py/.md, physical_derivatives.py/.md,
   detached_metric_wall_plasma.py/.md and WORKSPACE_MEMORY_DESIGN. The released
   pair still constructs its full velocity-gradient tensor before constructing
   and evaluating the full stress tensor. These independent evaluations retain
   both9-component tensors together. In a NEW opt-in subclass, move gradient
   construction until after stress_rate completes and stress is released. Keep
   each arithmetic expression, state, density, returned term, face contribution,
   physical mass and original check. This changes call scheduling, so do not
   reuse the earlier del-only AST-identity claim. Verify full coarse/fine old/new
   physical workflows and exact replays, and measure actual peak storage.
2. If rescheduling is insufficient, evaluate stress/gradient contractions by
   component or stream exact quadrature batches in another frozen-source variant.
   Preserve all cylindrical connections and force/work adjoints. Any changed
   floating-point reductions require actual physical component/work comparisons,
   not an assumption of equivalence or a relaxed bound. No floors or repairs.
3. Budget axial320+ and independent radial refinement. Radial field cells must
   subdivide the fixed64 fluid elements, so320 is the next admissible radial
   increment; it also increases the dense uniform metric more than axial
   refinement. Retain the improved self-terminating guard pattern and one large
   solver process at a time. The completed256→288 axial step does not replace
   doubled grids. Old B256/512 resource-interrupted attempt stays untouched.
4. Require full signed cases, independent q16/20 physical norms, each current
   component and actual Bt strong curl, face/local budgets and fresh full replay.
   Add amplitude and residual controls on new grids before declaring convergence;
   epsilon1e-5 sensitivity/raw-curl cancellation is not solved by these memory
   changes. Further fixed-field fluid refinement may be needed; current results
   do not separate every source of remaining Jz/radial-fluid error.
5. Continue full signed trajectory mesh/time/residual convergence and local
   budgets, then heated moving exact vacuum, transport/ionization/validity,
   actual prescribed coils and physical ends, converged formation and the UI.
   The radial1F8/64,B32/64 full100us control is still pending; near-grid growth
   in old trajectories is unresolved. Do not adopt new rate-only evidence as
   time stability, confinement or a laboratory experiment recommendation.

## Practical continuation

All current new files are frozen after calculation. Extend them rather than
editing a used Python source. Complete physical source identity for the detached
workflow is f2402996f3cdf2aca408a538ffaac8e850fea0ba964d7fb6ac8d1f385b1c43ae.
Baseline outputs: detached_current_fine; refined: detached_axial288; independent
current audits: detached_baseline_curl and detached_axial288_curl, all under
outputs/. Existing guarded_detached_refinement fixes F64/256 and accepts field
--nr/--nz; --original verifies a fresh same-source replay, not an algorithm change.
Use streamed_audit_comparison for whole-workflow bitwise old/new comparisons.

Two new manifests share496702418 nominal duplicate bytes from completed,
verified, byte-identical archived NPZs. Preserve all paths/hashes and never
modify either shared path in place. Check actual disk headroom before new runs;
no system-wide storage cause is asserted. No unrelated file was deleted.
Preserve outputs/frc_showcase/poster.png (user hash a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca).
Use bd --no-daemon and local sync --flush-only. No remote/publication. Provisional
lab envelope remains~1m quartz tube,10–20cm diameter,selectable gas,0.1–100mTorr
neutral fill and adjustable Helmholtz-like bias. Hardware stays prescribed
waveforms; none of these fixtures establishes3-D stability or confinement.

---

# Continue from verified tensor mass inverse and measured resource limits

All workers are terminal and collected; no tool handle needs resuming. The full
original FRC goal and all 12 acceptance beads remain open. No unit tests added
or run: continue with complete end-to-end physical workflows, independent
comparisons, convergence, conservation and full replays.

Read the newest GOAL_STATUS, FIELD_COMPONENT_REFINEMENT_DESIGN/VERIFICATION,
and outputs/field_component_diagnostics/checkpoint.json. The preceding09b9363
checkpoint and this continuation both made verified progress.

## Current evidence

- Full original/new fine baseline and component comparisons pass with fresh
  complete replays. Tensor inverse changes only factorization evaluation;
  original equations, masses, quadrature, boundaries and checks remain.
  Default geometry/time integrators remain unchanged. Coarse full workflows
  capture identical actual assembled operators; fine physical inputs agree.
- Radial current error persists at 5.96939%, despite tiny vector error.
  Fine original/tensor Jr difference is 1.04097e-5 of its reference norm.
  Worst strict residual 6.4967e-11 under1e-10; power defect <=8.75e-9 of
  exchanged powers. This is an instantaneous verification fixture, not formation.
- Existing fluid32/64 same-field comparison now resolves each current component:
  Jr errors11.6635%/5.96939%, Jz36.1204%/1.18700%. Full audit replays bitwise.
- Original B256/512 initialization ended with resource interruption, no signed
  input/rate. Worker36206 exit143 was collected. Never resume/restart from its
  partial summary alone or relabel that preserved run as a numerical failure.
- Tensor full fine baseline/replay peaks7.420/7.150GiB; observation replay7.432.
  CurrentTransport.pair first-sign cumulative peak4.829→6.558GiB; second signed
  electromagnetic pair7.061→7.432GiB. Nested high-water intervals are not
  exclusive/live storage. See tensor_resource_observation/calls.json and full
  verified replay. All successful tensor guards remained untriggered.
- Storage estimates3.094GiB forB256/512 and5.098GiB forB512/256 exclude large
  factor/workspace peaks. Do not infer doubled-grid feasibility from them.
- All37 new archives,318 preceding archives and22 older closures are verified;
  147-file source SHA8ae692975fec4c2fac1345f723b13bb428c1bb1affa2c72557f527d476c8b332.
  All historical calculation Python files and the user's poster are unchanged.

## Next actions: frc-7mm.6 and frc-7mm.7

1. Reduce actual temporary quadrature storage in a new opt-in module. Inspect
   current_transport.py and regular_current.py companions and the source before
   modifying anything. CurrentTransport.pair retains full velocity-gradient and
   stress tensors beyond their last use; advection also forms full gradient
   arrays. The second signed electromagnetic evaluation sets the observed
   global peak, so transport alone is not established as sufficient. Profile
   before/after using own-worker resources and unchanged physical checks.
   Do not edit frozen calculation modules; add a separate variant.
2. Require complete original/new signed coarse and fine workflows, identical
   physical inputs/operators where applicable, independent physical components,
   work, q16/20 and fresh full replays. Arithmetic changes need measured actual
   differences, not assumed bitwise equivalence. Keep all original failures,
   physical masses/density, basis modes, bounds and numerical budgets.
3. Retry fixed-reference F64/256 separate axial B256/512 then radial B512/256
   only with justified memory/disk headroom. Use the new self-terminating guard
   pattern; the old SIGINT-only guard failed to promptly stop native work.
   Require complete both-sign rates, component/strong-curl comparisons and full
   fresh replays. Never substitute total-current norm for Jr accuracy.
4. Continue full signed trajectory mesh/time/residual convergence and local
   budgets, then heated moving exact vacuum, transport/ionization, real prescribed
   coils/ends, converged formation and experimental UI. The radial1F8/64,B32/64
   full100us control remains pending; old near-grid growth is unresolved.

Shared-storage manifests cover only byte-identical completed archived replay
NPZs. Both paths now share immutable bytes: never modify either in place. No
research evidence was deleted. Check actual free disk before new allocations;
host-wide changes are not attributed to a particular process without evidence.
Preserve outputs/frc_showcase/poster.png (user hash a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca).
Use bd --no-daemon and local sync --flush-only. No remote/publication. Provisional
lab envelope remains~1m tube,10–20cm diameter, selectable gas,0.1–100mTorr neutral
fill and adjustable Helmholtz-like bias; hardware remains prescribed waveforms.
These fixtures establish neither 3-D stability nor confinement.

---

# Continue from separated curl arithmetic and current-component error

Read the newest GOAL_STATUS and CURL_EVALUATION_DESIGN/VERIFICATION.
**All workers are terminal and collected. No job needs resuming.**
Use complete workflows, independent physical comparisons, convergence,
conservation and full replays. Do not create or run unit tests. The complete
original FRC sandbox goal and all 12 acceptance beads remain open.

## Verified current evidence

- Original strict-rejection vectors reproduce bitwise. Simple projected-curl
  factoring is not a repair. Exact raw products/sums plus addition-low accounting
  leave 2.97e-15/RHS distributivity error, but actual dual Faraday remains
  5.62e-11 against 1e-11 on that saved iterate. No solver change is adopted.
- Full physical component spectra at epsilon 1e-3/1e-4/1e-5 reproduce original
  fields and replay completely. Jr errors 5.9751%/5.9694%/6.6721%; high-band
  norm errors 0.02750%/0.27966%/2.78699% approximately scale as 1/epsilon.
  The dominant |k|20 spatial-error hump persists. Above |k|64 contains 86.672%
  of the 1e-4 versus 1e-5 amplitude-difference energy, but only 0.2195% of
  the 1e-4 squared error. No accepted state was filtered. These are cell modes,
  not unaliased continuous harmonics. Full original integrals/Parseval verified.
- Differentiating the actual centered magnetic rate independently gives 5.9425%
  radial-current error versus 5.9694% in Jt, with 0.1108% weak/strong difference.
  This is an inside-cell volume comparison excluding distributional face jumps.
  Signed cross terms close the error identity to 2.16e-16. Both amplitudes replay.
- New complete signed epsilon 1e-3 F64/256,B256/256 rates pass every original
  check, max 1.91046e-11 under 1e-10. Full q16/20 physical audits and exact fresh
  replay pass. Increasing amplitude slightly worsens velocity/toroidal-field
  errors: numerical-noise reduction is not an accuracy or formation gate.

## Next actions: frc-7mm.6 and frc-7mm.7

1. Budget independent field-axial and field-radial refinement at the fixed
   radial16,N128 reference and fixed fluid representation. Measure each current
   component, actual Bt derivatives and local budgets. Do not adopt the tiny
   total-current norm as evidence for Jr accuracy. Use the epsilon 1e-4 baseline
   alongside amplitude controls; 1e-3 visibly changes nonlinear responses.
   Finest-mesh memory and sparse-factor storage need explicit accounting.
2. Investigate accurate raw-curl and electric-sum representation. Applying C
   then P alone did not solve the original rejection. Keep original physical
   masses/density, real equations, adjoints/work, kernel terms, all residuals
   and bounds. A proposed arithmetic implementation needs complete exact-input,
   physical-power/field comparisons and signed full replays before adoption.
3. Continue complete signed radial1 F8/64,B32/64 100us trajectory controls
   against preserved F4/64,B16/64 when accuracy gates justify adoption, including
   time/mesh/residual and independent local conservation. Near-grid growth is
   unresolved; do not treat new instantaneous rate audits as time convergence.
4. Continue heated moving exact vacuum, transport/ionization and validity,
   physical prescribed coils/ends, converged formation and experimental UI.
   The user's ~1 m tube, 10-20cm diameter, variable gas, 0.1-100mTorr neutral fill
   and adjustable bias remain provisional. Hardware stays abstract waveforms.

## Evidence and preservation

- outputs/curl_evaluation_probe and its _replay/_verification: all four original
  iterates, alternative full actions and final exact raw reference.
- outputs/component_axial_spectrum and its _larger variant, both with complete
  _replay/_verification: full physical component spectra and amplitude effects.
- outputs/current_rate_curl_comparison and _replay/_verification: independent
  magnetic/current-rate decomposition, q16/20 and actual field fingerprints.
- outputs/radial16_amplitude_larger and _replay/_verification: new full signed
  physical rate workflow. Every coefficient/face array and metric replays.
- outputs/curl_spectrum_report_layout/verification.png: visually inspected.
  Original report layout also preserved. No PDF generated in this continuation.
- outputs/curl_evaluation_diagnostics/checkpoint.json: complete source, input,
  worker and replay inventory. No disk guard triggered. Large geometry builds
  and replays used separate processes.

Five verified groups comprise four diagnostics and one full signed rate
calculation. Six new Python modules are frozen, each with its companion.
18 new archives, 300 earlier archives and 22 earlier closures verified.
Final 140-file source SHA b0ee5f17e5baf0fa4fc564d3fe40437173b05479d16210e2ef040f0a1d0ce922.
No old calculation source changed. Preserve user poster SHA
 a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca
outside research commits. All work remains local; no remote/publication.
Older handoffs below are historical.

---

# Continue from component accuracy and split-residual path diagnosis

Read the newest GOAL_STATUS and SPLIT_RESIDUAL_PATH_DESIGN/VERIFICATION.
**All workers are terminal and collected; no job needs resuming.**
No unit tests: use complete workflows, independent physical references,
convergence/conservation measurements and full replays. Keep all 12 acceptance
beads and the overarching FRC goal open.

## Current evidence

- The original finest tighter-bound failure is reproduced with unchanged input,
  equations, residuals and budgets. New observations show final opposing-vector
  gap 6.87234e-11 of RHS. Electric curl distributivity leaves only 0.4245% of that
  gap; H/T and current/inertia paths do not dominate. All-term remainder is
  5.41e-17 of RHS. The complete fresh standalone observation replay is bitwise.
- Exact Decimal sums verify 786432 actual electric coefficients, with zero
  mismatches and a complete postprocessing replay. Propagating addition loss
  alone leaves 66.21% of the gap unexplained. Do not conflate addition rounding
  with the full curl-action discrepancy or infer an accepted solver repair.
- Complete independent physical component audit/replay on F64/256, B256/256,
  fixed radial16, N128 reproduces every original field/reference fingerprint.
  Jr-rate errors are 0.05969389 / 0.06672077 at epsilon 1e-4 / 1e-5. Jz errors
  are 0.01187001 / 0.01193623; Jtheta error is 0.0001904658. The total J error
  0.000878993 therefore does not certify these smaller physical components.
- Radial-current amplitude difference is 0.03014823 of its own continuum norm.
  88.674% of its squared difference is interior and 11.115% in the outer 5% R.
  Radial-current error itself is about 6% in each band. Axial-current and radial-
  velocity squared errors are more wall-localized (58.08% and 72.05% there).
  All nine fluid/magnetic and all E/J components are measured; q16/20 norm and
  masked-band differences are recorded. Band integrals mask original nodes,
  rather than subdividing cells at the two fixed physical band edges.

## Next actions: frc-7mm.6 and frc-7mm.7

1. Inspect original `StrongWallGeometry.faraday = projection @ original.faraday`.
   On the exact captured original iterates, compare assembled action with raw
   curl followed by the same projection, and isolate product/sum rounding.
   Keep actual electric terms and the independently verified addition low part.
   Do not replace a rejected full rate with a diagnostic vector. Any proposed
   arithmetic change must independently verify the original full equations,
   physical fields and work, preserve strict checks/budgets, and fully replay
   complete signed workflows. The observations are not a repair.
2. Determine whether component errors arise from representation, axial/radial
   differentiation, nonlinear amplitude effects or arithmetic. Use physical
   component norms and fixed independent radial16, N128 reference. Axial-mode
   diagnostics or separate directional refinement can guide the next grid;
   do not infer a cause merely from the dominant total-current norm. A causal
   link to the strict residual-path mismatch remains unproven. Any Fourier
   analysis is diagnostic only and must not filter accepted states.
3. Continue complete signed radial1 trajectories at F8/64, B32/64 against the
   preserved F4/64, B16/64 runs when the accuracy gates justify adoption. Require
   100 us, time/mesh/residual controls, local budgets and independent physical
   conservation. The prior near-grid growing branch remains unresolved.
4. Continue heated moving exact vacuum, defensible transport/ionization and
   validity diagnostics, real prescribed coils/ends, converged formation and
   the experiment UI. Keep original lab inputs provisional and distinguish
   reversal, closed flux, equilibrium, 3-D stability and confinement.

## Artifacts and resource constraints

- `outputs/split_residual_probe`: original full exact-input observation.
- `outputs/split_residual_probe_replay_v2` and `_verification`: passing complete
  standalone replay. The older `_replay` directory is partial and preserved.
- `outputs/split_residual_arithmetic_audit` and component_amplitude_localization:
  original audits plus their `_replay` and `_verification` directories.
- `outputs/split_precision_report/verification.png`: visually inspected
  scientific figure, with a PDF alongside. PDF not separately inspected.
- `outputs/split_residual_diagnostics/checkpoint.json`: source/input/worker
  inventory, resource guards and preserved interruption.

Run large original/replay geometry builds in separate processes. The earlier
in-process replay was interrupted when free disk fell below 256 MiB; no global
cause is asserted and no unrelated files were deleted. The standalone replay
completed with the same physical calculation and guard. Source-derived finest
trajectory modal storage still exceeds 15 GiB before field/work arrays.

All six new modules and old calculation sources are frozen; read companions
and extend via separate modules. Final 134-file source union:
150dc2f3e7614eed5da58d314b6ff3e8cf804b641c1f9fec9f91b738e4712745.
Three complete replay groups and 12 new archives verified; 288 preceding
archives and 22 older source closures unchanged. User poster SHA
 a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca.
Preserve it outside research commits. Beads and research remain local; no
remote or publication. Older handoffs below are historical.

---

# Continue from complete joint refinement and physical amplitude sensitivity

Read the newest GOAL_STATUS, JOINT_E2E_DESIGN and JOINT_E2E_VERIFICATION.
User policy: no unit tests; use complete workflows, independent physical
comparisons, convergence/conservation measurements and full replays.
**All workers are terminal and collected. No job needs resuming.**
All 12 full-goal beads remain open, especially frc-7mm.6 and frc-7mm.7.

## Verified evidence

- Fixed radial16, N128 reference SHA
  3b142fe986aca23d340b1739c17ec83b5b94e04574cff2a661fe1b0a8527ef08.
  Full F32/256 and F64/256 signed cases on B256/256 and full reruns pass.
  Finest physical rate errors: vr 0.00183142646, Br 0.00097962351,
  Btheta 0.00547509477, Jt 0.00087899301, E 0.00000456807405.
- Finest epsilon 1e-5 success and complete replay give Jt 0.00088543967.
  The centered current coefficient difference is 0.00947096, while actual
  physical vector difference is 0.000107450107. Radial/axial current-component
  physical differences are 0.0301482272 / 0.00126286432 relative to their own
  references. E difference is 7.148893e-8. The independent reconstruction
  reproduces original physical fingerprints and fully replays at q16/20.
- The tighter 1e-11 case rejects the identical minus input and reproduces its
  rejection bitwise. No failed rate exists. Saved raw input:
  outputs/radial16_finest_controls/case_01.npz, SHA
  ce34210139b000f4c2746dfde664ba1dd188d5e2485572a364e184bf12c4d999.
  Dual Faraday residual stays near 5.89e-11; scalar norms imply residual-path
  mismatch at least 2.499e-11 at the recorded final iterates. This does not
  identify the cause or prove impossibility at other iterates.
- Both separate radial1, N64 signed 100 us / 60-step trajectories, every saved
  stage/state and all physical/conservation audits replay exactly. Actual
  energy drift is at most 3.636e-7 of initial physical signal-energy scale.
  These coarse integrations do not establish refined time/mesh convergence
  or remove the previously observed near-grid growing branch.
- 20 complete replay groups (including postprocessing, rejected cases and
  completed-case selection), 137 raw-file comparisons and 71 source archives
  verified. The preceding 42/23/59/61/32 archives and 22 older source closures
  remain unchanged. See outputs/radial16_joint_diagnostics/checkpoint.json.

## Next concrete actions

1. Reconstruct the exact rejected input with its current split-inverse and
   direct-metric setup, enforcing bitwise input identity. Capture the actual
   operator/dual residual vectors and individual mass-inverse actions without
   changing the arithmetic or acceptance. An observation-only copy must
   reproduce the original first-attempt rejection exactly before interpreting
   its decomposition. The old balanced WallRoundoffCore is different arithmetic
   and cannot be used as a drop-in observation of this failure. Preserve every
   original bound, physical mass, actual density and refinement budget.
2. Investigate the radial-current amplitude sensitivity in physical component
   norms, including near-axis/wall structure, and reconstruct fluid/magnetic
   difference fields as needed. The existing current/E audit measures these
   fields directly; close total error norms alone do not certify sensitivity.
3. Continue full signed radial1 trajectories on F8/64, B32/64 against preserved
   F4/64, B16/64 results, keeping the independent reference and physical inputs
   fixed. Require full 100 us runs, physical field/conservation audits, complete
   replays and time controls. Preserve failures and growth. This candidate grid
   is not presumed accurate. Budget memory first: source-derived retained dense
   modal storage is 0.0579 GiB here, but 15.414 GiB on F64/256, B256/256, before
   field matrices, sparse factors, nonlinear/Krylov and audit workspaces.
4. After those gates, continue local current/enthalpy/Poynting budgets, heated
   moving exact vacuum, transport/ionization, actual coil/end conditions,
   converged formation and the experiment UI. Keep the full FRC goal open.

## Preservation and local state

All calculation code used here is now frozen, including the 14 new opt-in
modules listed in source_preservation.json. Read companions before extending
it; use new modules rather than editing archived sources. The final 128-file
source union is archived under radial16_joint_diagnostics/final_source.
Earlier interrupted batches and partial replays stay intact; completed-case
selection explicitly preserves the incomplete parent status. The finest
control disk guard did not trigger. No unit tests were added or run.

Scientific figure: outputs/joint_e2e_report/verification.png, visually inspected;
PDF is alongside it. User poster SHA
 a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca.
Do not include that preexisting poster change in research commits. Work and
beads are local with no configured remote or publication. The provisional lab
envelope is not a selected experiment; no 3-D stability or confinement claim.
Older entries below are historical.

---

# Continue from directional term accuracy and strict arithmetic recovery

The full FRC experiment-sandbox goal and all12 acceptance beads remain open.
**470 tests passed in233.11s. Every worker is terminal and collected.** Read the
newest GOAL_STATUS.md, DIRECTIONAL_TERM_DESIGN.md and
DIRECTIONAL_TERM_VERIFICATION.md before changing code. No local or remote job
needs continuation or restart.

## Current evidence

- Same physical radial16,N128 reference, growth12307.1561445/s; original reference
  array SHA3b142fe986aca23d340b1739c17ec83b5b94e04574cff2a661fe1b0a8527ef08.
  Axial doubling improves Br .051938 to .00231942 but leaves J .0847343.
- Radial field128->256 at fluid32/128 lowers vr .1088365->.0309084 and
  J .0875688->.0337504. Br .051625 and Btheta .279381 remain large.
  Fluid32->64 alone improves pressure error31.53-fold but leaves magnetic-force
  error. Refine the two approximations independently rather than attributing
  the full error to one term from a single net residual.
- Independent original pressure/Lorentz and ideal/Hall/pressure/inertia terms:
  full captured fluid sum bitwise matches saved rates; physical terms/curls
  checked independently. Hall/pressure toroidal norms ~24000 times net response.
  Inertia includes numerical-current feedback. Raw errors, cross integrals,
  strong/projected curls, fixed physical wall bands and q16/20 norms remain.
- Original radial256 solve failed with operator2.60e-10 and dual2.83e-10 under
  1e-10. Original failure and exact input77e43521... remain preserved. HT
  association was ruled out; mass-inverse distributivity dominates opposing
  residuals. Independent reconstruction retains7.06% vector remainder.
- Opt-in SplitInverseWallPlasma/Core uses E=M^-1f+M^-1(IJt), sharing forcing
  inverse with original RHS. Every original equation, physical mass, actual
  density, boundary, operator, bound and budget remains. Full exact input passes
  at5.27777e-11, bitwise replayed. Polynomial full B/J/E/work and atomic tests pass.
  Full power defect3.99063e-9 and restored full Faraday rate-relative2.722e-7
  remain explicit; solving accepted equations does not certify physical accuracy.
- Three full signed split-inverse pairs replay bitwise: baseline, radial256
  eps1e-4, radial256 eps1e-5. Physical errors nearly agree across amplitudes;
  centered current coefficient difference .00108516 is retained. The old/new
  baseline inputs and full fluid rates match bitwise; J coefficient change
  5.86e-10. The guarded trajectory sourcef255a576 remains unchanged/unadopted.
- Thirteen complete replay groups,42 new source archives verified;23/59/61/32
  preceding archives byte-unchanged and22 earlier source closures unchanged.
  Full checkpoint: outputs/radial16_directional_diagnostics/checkpoint.json.
  Visually inspected report: outputs/radial16_directional_report_v2 PNG/PDF.

## Immediate next action: frc-7mm.7 and frc-7mm.6

1. Read companions for split_inverse_wall_core/induction/plasma,
   split_inverse_radial_rates, radial_parent_terms, radial_term_capture/audit,
   radial_drive_capture/audit and radial_error_localization. All now have frozen
   source closures; do not edit code imported by archived or live calculations.
2. Run a fresh combined-grid comparison F32/256,B256/256 and F64/256,B256/256,
   retaining fixed independent radial16,N128 reference and physical masses.
   Use both signed amplitudes1e-4 and1e-5 where possible, q16/20, strict original
   residuals and explicit failed records; full fresh replay is required. Separate
   initial representation, complete rates, pressure/Lorentz cancellation,
   current and one-sided face errors. Never make a rejected state look accepted.
3. Localize any remaining pressure/current/wall error before adopting the new
   arithmetic into complete growing-response trajectories. Then repeat actual
   signed100us mesh/time/residual convergence and continuous conservation/local
   enthalpy/Poynting/charge face budgets. Old time order remains unestablished.
4. Continue heated moving exact vacuum, transport/ionization and validity,
   physical ends/coils, converged formation and experimental controls/UI.
   No atmosphere, filtering, artificial mass, heat or energy repair.

The roughly1m quartz tube,10-20cm diameter, selectable gas,0.1-100mTorr neutral
fill and adjustable Helmholtz-like bias remain provisional lab variables.
The current smaller argon filled-periodic case is a numerical verification
fixture, not an experiment recommendation. Axisymmetric closed flux does not
establish three-dimensional stability or particle confinement.

No remote is configured; all research and issue state stay local. Preserve
outputs/frc_showcase/poster.png (SHAa783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca).

---

# Continue from fixed radial16 accuracy and strict solve recovery

The complete FRC experiment-sandbox goal remains open. Read GOAL_STATUS.md,
RADIAL_MODE_DESIGN.md and RADIAL_MODE_VERIFICATION.md first. **462 tests pass
in 230.35 s.** The active new entry point is examples.balanced_radial_rates.run
for further instantaneous studies; do not silently change guarded trajectories.
Older handoff sections below are historical, including their old next actions.

All jobs from this continuation are terminal and collected. Seven full bitwise
verification groups and23 new source archives are verified. The59 preceding,
61 older matched and32 older wall archives are byte-unchanged, and22 earlier
verified source closures still match. No worker needs resuming. The full
checkpoint is outputs/radial16_diagnostics/checkpoint.json.

## Current scientific evidence

- Fixed independent radial16 reference at N128, checked against N256. Growth
  12307.15614451037/s; full physical/J/E component shape differences <6.9e-9.
  Raw reference outputs/radial16_reference/continuum_01.npz SHA256
  3b142fe986aca23d340b1739c17ec83b5b94e04574cff2a661fe1b0a8527ef08.
- The original N64 archive remains bitwise unchanged but has a .003251 off-grid
  curl/current mismatch from its real-Nyquist convention. Do not use it as a
  resolved high-radial reference. New N128 independent curl checks pass.
- Separate fluid/field/axial meshes show large full-RHS refinement effects.
  F16/128,B64/128 vr/Br/Jt errors3.86912/.129238/5.10160 improve to
  .1088365/.0519380/.0875688 at F32/128,B128/128. **Btheta error .383069**
  remains; this is not an accepted spatial accuracy or formation gate.
  Wall-region vr/Btheta/Jt errors .576775/1.154478/.438075 identify a concrete
  next diagnostic target; these normalize to the local reference norm.
- Initial final n/vr errors~1.05e-5, Bz .000116518, weak J .000266164;
  radial normal-J/tangential-B jumps .000320212/.000323249. Axial jumps
  remain9.93422e-5/.000126339. Strong Ampere .000982571 uses strong J norm,
  excludes distributions and differs from the old extracted-mode denominator.
- The original F32/B128 minus solve strictly rejects at operator1.39069e-10
  versus1e-10. Its raw input SHA256 is
  e7a3cf637ceff0f3201a88f12ddd84b593f28d23f898406915c8d040a3eac755.
  BalancedWallCore proposes the mean of operator/physical dual corrections
  when the operator bound fails. Both original residual paths then pass;
  worst7.72779e-11. All physics, tolerances and refinement budgets unchanged.
  Exact saved-input recovery and baseline bitwise equivalence are recorded.
- Independent analytic cylindrical B/J/E and work tests pass on6/12/24 cells;
  iteration failure is atomic. Full actual power cancellation on the saved
  input is2.696e-9 relative to absolute terms, separately measured.
- Amplitude1e-4/1e-5 gives nearly the same physical error norms; centered
  coefficient Jt differs .000773103 relatively. Keep that limitation explicit.
  Volume quadrature16/20 norms differ <=3.05e-11 relatively. The largest
  axis relative discrepancy belongs to7.84e-11 error and is1.12e-15 absolute.
- Inspected scientific PNG/PDF: outputs/radial16_report/verification.
  Full inventory: outputs/radial16_diagnostics/checkpoint.json.

## Concrete next actions

1. Read companions for archived_contact_mode, radial_contact_fixture,
   balanced_radial_rates, radial_mode_interfaces, balanced_wall_core,
   balanced_wall_induction and balanced_wall_plasma. These are opt-in; all
   earlier source closures and current accepted trajectories remain frozen.
2. Hold the same radial16,N128 reference fixed. Next compare axial refinement
   dict(elements=32,points=256,field_nr=128,field_nz=256,epsilon=1e-4)
   separately from radial-field refinement
   dict(elements=32,points=128,field_nr=256,field_nz=128,epsilon=1e-4), then
   refine fluid elements to64. Account for quadrature20 memory before launching
   concurrent large cases. Every destination must be new and immutable.
3. Localize remaining vr/Btheta/Jt rate error by original pressure/current/work
   contributions and physical wall/interface regions. Retain actual inputs,
   signed rates, full E/J, independent quadrature, amplitude, strong/weak/face
   errors and exact rejection. Do not fit fields or remove species feedback.
4. Once spatial accuracy supports it, repeat complete signed100us response
   mesh/time/residual controls with local budgets and small-signal energy scale.
   The old evolution source is stillf255a576 and has not adopted balanced
   induction. Previous medium final J error .0431611 and unestablished time
   order remain active limitations; rate tests do not supersede them.
5. Continue full heated moving exact vacuum, transport/ionization, real coil/end
   boundaries, prescribed-waveform formation and experiment UI. All12 unfinished
   beads, including frc-goal, remain open until actual acceptance evidence.

Use the repository .venv with PYTHONDONTWRITEBYTECODE=1,
VECLIB_MAXIMUM_THREADS=1, MPLCONFIGDIR=/private/tmp/plasma-bench-mpl,
MPLBACKEND=Agg. New rates use examples.balanced_radial_rates.run; repeat from
scratch and verify with examples.wall_diagnostic_replay.verify. Do not modify
any imported source of a running or frozen calculation. Use bd --no-daemon;
sync issue state locally. No remote is configured; no publication is authorized.

Preserve the user's modified outputs/frc_showcase/poster.png. The proposed lab
remains~1m quartz tube,10–20cm diameter,0.1–100mTorr neutral fill, selectable gas
and adjustable Helmholtz-like bias. No lab design,3-D stability or confinement
has been established. This remains useful ongoing research, not a blocked goal.

---

# Continue from the radial near-grid response failure

The complete FRC experiment-sandbox goal remains open. Read the newest
GOAL_STATUS.md, PHYSICAL_INTERFACE_AUDIT_DESIGN.md and
PHYSICAL_INTERFACE_VERIFICATION.md. **457 tests pass in254.80s.** No production
physics/time-solver code was edited. New source closures are opt-in diagnostics.

## All jobs terminal

Both larger-grid60-step/100us primary runs and complete reruns finished; every
accepted state/stage and non-runtime metric compares bitwise. Their verification
records are outputs/matched_mode_guarded_ne4_z64_t60_plus_verification and the
corresponding minus_verification. Former worker sessions10675/67242 exited0 and
were collected; do not restart or poll them. The coarse120-step pair, both full
reruns, every physical audit, history, interface, extracted-mode diagnostic and
regression are also terminal. No job was killed or accepted with relaxed bounds.
All59 new source archives are verified;61 preceding and32 older archives plus22
earlier verified source closures remain unchanged. See the full checkpoint in
outputs/physical_interface_diagnostics/checkpoint.json.

## New evidence

- Four additional complete100us signed runs, coarse120 steps and medium60steps.
  Medium F4/64,B16/64 intended projection improves to2.40217, but final J error
  is.0431611 and Bz.00496073. This is a failed accuracy gate; no mesh convergence.
- Coarse30/60/120 controls do not establish fourth-order time convergence.
  Representative uniform spectra have omega*dt118/318 at60steps on coarse/
  medium grids. Phase underresolution is a risk, not the sole diagnosed cause.
- Full physical histories replay every saved checkpoint atq12/16, matching the
  separate finalq16 audit. Medium extra current grows~71682/s after cancellation
  of the independent slow reference; direct actual spatial RHS projects~71426/s.
  Probe amplitude1e-4/1e-5 current rate difference1.57e-6; q16/20 agrees.
- Extracted current is97.83% azimuthal, B98.07% axial; axial Fourier1 dominates.
  Radial Bessel15/16 dominate density/Bz on16 field cells.32 modes capture~94%
  of density/Bz. Weak J versus strong volume curl(B)/mu0 differs.599823; this
  excludes distributional face terms. Independent matched continuum action has
  large radial/induction discrepancies; no validated physical eigenmode.
- Exact fluid subdivision to8/16 elements preserves fields<7e-14 and leaves
  current rate nearly unchanged, but changes momentum rate dramatically. It is
  a separate diagnostic, not an established repair. All three extracted-mode
  studies fully replay every array and non-runtime metric.
- Eight actual original physical face/volume budget cases and axis/wall audits
  replay bitwise atq16/20. Energy cell L2 defect.60485 ->.00017354; final
  particle/current/eentropy/ientropy defects1.306e-7/3.791e-7/9.741e-5/4.022e-8.
  Energy absolute trace jumps5.19e-4 ->1.84e-7 W. Max quadrature norm discrepancy
  9.64e-10. These are rate consistency results, not full nonlinear acceptance.
- Initial conservative energy Hessian is independently verified against70digit
  original EOS differentiation and analytic derivatives. Actual total-energy
  drift/initial perturbation scale<=6.42e-7 across five signed pairs; the scale
  is not an independently conserved free energy. All energy audits replay.
- Inspected report: outputs/physical_interface_report/verification.png and PDF.
  New full regression: outputs/physical_interface_diagnostics/pytest_final.log.

## Required next work

1. Resolve the near-grid radial branch using independently matched fluid/field
   radial spaces and fixed-physical-wavenumber refinement, retaining exact
   physical masses, full species/current/work coupling and all wall constraints.
   Concrete first control: the independent radial16,N64 mode is already frozen
   in outputs/stationary_contact_audit/continuum_07.npz (SHA256
   f9c92d8d94be18b10955fdde1f20868e2d3f0fafcc58a7157eca3e318258f6c6),
   with growth12307.156241/s. Add an opt-in reference/fixture that accepts this
   archived mode without changing continuum_contact_mode or matched_mode_fixture.
   Compare actual full rates on field radial16/32/64 (axial64 then128), varying
   fluid elements4/8/16 separately and keeping the same physical radial16 mode.
   Preserve exact initial representation, full E/J, strong and face errors and
   epsilon1e-4/1e-5 controls; refine the independent axial reference separately.
   A comparison that simply follows a different largest eigenvalue as the grid
   changes is insufficient. Include strong/weak/face consistency and uncaptured
   norms; no fitted eigenmode, filter, artificial mass or energy repair.
2. Repeat complete growing signed trajectories with genuine mesh/time/amplitude/
   residual convergence. Keep small-signal energy and actual local interfaces.
   Source for current guarded evolution remainsf255a576; extracted diagnostics
   do not change the accepted trajectory solver or physical model.
3. Continue heated moving exact vacuum (never an atmosphere), transport and
   ionization validity, real coil/end boundaries, prescribed-waveform formation,
   experiment controls/UI. Keep frc-goal and all current acceptance beads open.

Preserve outputs/frc_showcase/poster.png, which was already user-modified. No
remote is configured and no research data was published. The user's~1m quartz
tube,10–20cm diameter,0.1–100mTorr neutral fill, selectable gas and adjustable
Helmholtz-like bias remain exploration inputs. No laboratory formation,
three-dimensional stability or particle-confinement claim follows.

---

# Continue from matched growing-response verification

The full user goal remains active. Read the newest GOAL_STATUS.md,
MATCHED_MODE_SOLVER.md and MATCHED_MODE_VERIFICATION.md. **446 tests pass in
234.86 s.** All new jobs are terminal. All six signed trajectories replay every
array and non-runtime metric bitwise. See outputs/matched_mode_diagnostics/checkpoint.json
for archive and final verification details. No running job needs continuation.
The overall scope has not changed.

## New evidence

- Eight independent field-direction rate cases, all arrays/physical fingerprints
  and both quadratures bitwise replayed. Final Br/Jt errors0.000244246/0.000219551
  at field64/128, fluid16/128; source 7bb77514. The old0.55/0.48 gate is resolved
  as an instantaneous accuracy study, not full fine-grid trajectory validation.
- Modal uniform physical block inverse and physical uniform current-inertia
  preconditioner preserve all actual equations, full Fourier modes and strict
  bounds. Focused independent operator/stage/atomic tests and full regression.
- Six full 100 us signed runs on coarse fluid 4/32 and field 8/32: 30/60 steps at
  epsilon 1e-4 plus30 steps at2e-4. Source f255a576. Three complete independent
  physical audits and two physical time/amplitude comparisons replay bitwise.
- At60 steps n/Br/J/E errors 0.001011/0.001937/0.002748/0.008678. Growth projections
  near 2.402 agree with independent 2.402175837, but initial representation is
  explicit and this alone is not convergence or a confinement result.
- Half-dt changes E by 0.00136156; double amplitude changes all measured responses
  <=2.07e-6. Energy normalized by total background is near roundoff, not a
  small-signal conservation proof. Full physical quadrature 16/20 norms agree.
- GuardedModalGauss bounds unproductive proposal work without changing acceptance.
  Prior slower pilot algorithms and one failed negative trajectory are preserved.
  Eleven pilot accepted stages/endpoints have fresh independently replayed audits.
  A stale-import postprocessor failed before producing physical audit output;
  its frozen on-disk source limitation is explicit in its termination.json.
- Report PNG/PDF visually inspected. All 22 old verified calculation closures and
  all 32 preceding checkpoint archives unchanged. Preserve user's modified
  outputs/frc_showcase/poster.png. No remote or publication.

## Immediate next action — frc-7mm.7 remains open

1. Complete full response mesh/time refinement. Start with
   examples.guarded_mode_evolution.run on spec
   dict(elements=4,points=64,field_nr=16,field_nz=64), both epsilon +/-1e-4,
   duration1e-4 and60 steps. Then independent matched_mode_audit and full
   matched_mode_replay. Every step has an immutable accepted checkpoint.
   Finer target: field32/128, matching fluid axial128, with radial element
   refinement. Fine rate results do not replace actual nonlinear trajectories.
2. Add a 120-step coarse pair to quantify temporal convergence; use the physical
   comparison, not merely a fitted amplitude. Test tighter residuals separately
   and preserve any rejection. The original extra-tight1e-14 prototype rejected;
   never silently relax bounds or claim it passed. Measure initial small-signal
   energy scaling before interpreting background-normalized roundoff drifts.
3. Close actual local current/enthalpy/Poynting interface and jump budgets.
   Continue heated moving exact vacuum (no numerical atmosphere), transport/
   ionization and model validity, physical ends/coils, converged formation and
   experimental controls/UI. Keep all acceptance beads and frc-goal open.

All analysis is fully ionized bounded periodic verification. The user's roughly
1m quartz tube, 10–20cm diameter, 0.1–100mTorr neutral fill, selectable gas and
Helmholtz-like adjustable bias are provisional exploration variables. Closed
flux or axisymmetric response does not establish 3-D stability/confinement.

---

# Continue from the wall-space correction and matched continuum gate

The full experiment-sandbox goal remains active. Read the newest GOAL_STATUS.md,
STRONG_WALL_DESIGN.md and STRONG_WALL_VERIFICATION.md. **436 tests passed in
232.03 s.** All new simulation, full replay and regression-test jobs are terminal.
No local or remote run needs continuation or restart.

## Current evidence

- Original rapid mode is localized near the outer wall with incompatible
  normal flow/Btheta/current traces. Independent physical reconstruction and
  radial/Fourier projections retain uncaptured norm and coarse background error.
- New explicit wall Galerkin curl Fw=P F, inactive wall-normal fluid momentum,
  admissible magnetic Krylov vectors and nonlinear Gauss stage coordinates.
  Strict new equations, full old raw-curl defect and zero reaction work are
  explicit. Original source classes and archived failures are preserved.
- Actual 344-coordinate contact growth 9203.373408/s matches the separately
  restricted original Jacobian. Three full matrices and all physical/projection
  arrays replay bitwise. Source ded3bac43979f5bde181c2eb6989d6af28b57f1d64bca09dd3f94413b7203a9a.
- Five new complete 100-microsecond contacts and six full waves; source
  d7a422b0f3e147e84e1acbb5b3db5c17b4eb162eab5d2d55845d829e18ebbc96.
  Wave time ratios 15.967079/15.991774; contact 15.883252. Contact energy drift
  decreases fourth order. Finest contact n error 9.97995e-7, pressures ~3.8e-6,
  E 5.00415e-4; unwanted B 1.00896e-13 T retained. All arrays/non-time metrics
  replay bitwise, with quadrature 16/20 audits.
- Independent original-continuum growing mode rate comparisons at mesh 1/2/4/8/16,
  quadrature 16/20 and amplitude 1e-4/1e-5: all 12 cases replay bitwise. Source
  c100c3eea05acbb767c65887bf692335ff88b788d884eef3839aaf7e5b0ec4fb.
  Finest full E rate comparison error 0.00648538, radial velocity 0.002272;
  radial-B/current rate errors remain 0.547444/0.480471. Do not claim adoption.
- 22 prior verified source closures and 32 previous-checkpoint archives unchanged.
  All 32 new source archives verified. Inspected scientific report PNG/PDF:
  outputs/strong_wall_report; checkpoint: outputs/strong_wall_diagnostics.
  Preserve user's modified outputs/frc_showcase/poster.png; no remote/publication.

## Immediate next action — frc-7mm.7 remains open

1. Read companions for continuum_contact_mode, continuum_wall_probe,
   strong_wall_geometry/plasma/induction, wall_induction_core and Gauss coordinates.
2. Refine radial and axial field meshes separately for the independent continuum
   mode. Current fixture couples default field resolution to fluid elements;
   generalize in a new explicitly frozen protocol, retaining analytic state,
   full physical rates/E/J, initial error and axis/wall/max/quadrature diagnostics.
   Slow magnetic response requires cancellation accuracy far beyond a B-field
   representation or fast-wave check. Do not suppress the residual mode.
3. Independently evolve small continuum perturbations with time/amplitude/mesh
   controls; compare physical fields/growth, not differing largest eigenvalues.
   Complete local charge/enthalpy/Poynting interface budgets including jumps.
4. Then continue heated moving exact vacuum, transport/ionization, actual coil/end
   boundaries, converged formation and experiment UI. Keep the complete goal and
   numerical acceptance beads open until their actual criteria pass.

This is substantial progress, not a blocked goal or a reduced objective. Older
sections below are historical; do not restart old terminated inaccurate jobs.

---

# Continue from implicit verification and exposed contact-growth gate

The full experiment-sandbox goal remains active, with twelve unfinished beads.
Read the newest GOAL_STATUS.md section, IMPLICIT_PLASMA_NUMERICS.md and
IMPLICIT_PLASMA_VERIFICATION.md. **432 tests passed in 227.47 s.** All jobs from
this checkpoint are terminal, including deliberately stopped inaccurate
contact candidates. No local or remote job needs continuation/restart.

## Current evidence

- Midpoint/Gauss/Radau advance the actual full physical-mass plasma equations
  with unchanged original induction bounds, fixed physical coordinate scales,
  strict nonlinear equations and atomic rejection. A fixed-background block
  preconditioner plus optional matrix-free Newton solves the equations without
  substituting the background into evolution. Dense coupling storage is only
  suitable for bounded verification, not large formation grids.
- Eleven full waves replay all arrays/non-time metrics bitwise; independent
  quadrature16/20 includes full E, both temperatures and physical energy.
  Time ratios: midpoint3.975482; Gauss15.967079/15.991775;
  Radau7.975201/7.994196. Radau has declared measured numerical temporal damping.
  Source feeb358468dfd0074e61effcca9af8dcb727748f113907d22b4e0e74e958c8ef
  (Gauss/midpoint),445623fdadd9f0825ae5441dc4c3fe26e92859e5c5d53564d5639949c9bb4240
  (Radau). Report: outputs/implicit_plasma_report/convergence.png and PDF.
- Three long Gauss contacts and one Radau contact developed large errors and
  were stopped, with terminal records and saved accepted states. Do not claim
  an accurate full100us transit. Initial chord-only failed pilots remain frozen;
  the captured ninth-step Gauss failure now has a passing Newton regression.
- A 360-coordinate stationary-contact spatial Jacobian has real growth
  854532.963/s at epsilon1e-5 versus854532.845/s at1e-4. The uniform rest control
  has maximum real4.28e-9/s. All matrices replay bitwise. Full source
  b82a7bf386a9d591e67699db71c9a5326e35c4bf709a30364b208645edd2449f.
  The uniform null-mode relative mixed error is ill-conditioned; separate
  nonzero independent rest-operator comparisons are <=4.403e-12.
- Paired species-pressure/current-entropy feedback is needed for rapid growth
  in saved-matrix ablations. Removing it is a diagnostic, not a solver fix.
- A separate continuum Fourier/Bessel operator imports no production spatial
  solver, verifies three uniform wave branches and has slower growing contact
  modes. First-mode axial refinement converges to8763.749254/s; higher sampled
  Bessel modes give about12,000/s. Eight matrices replay bitwise. Audit source
  b7dde0bbad1b307e752bc2a9b441e158dd7c790dd88728c684782667aee647a2.
  **The radial spaces differ: no matched-mode or numerical/physical cause has
  been established.** Spectral growth alone is not a laboratory stability claim.

## Immediate next action — frc-7mm.7

1. Read companions for contact_linear_probe, stationary_contact_continuum,
   stationary_contact_audit, regular_species, species_projection, regular_plasma
   and the relevant independent physical sampling code. Saved production modes
   and continuum matrices contain everything needed for the next comparison.
2. Reconstruct the production leading eigenvector as physical density, velocity,
   both entropy/pressure fields, B/J and local charge/face jumps. Compare its
   matched radial/axial subspace against the independently linearized parent
   equations. Consider an explicitly divergence-free continuum magnetic basis
   to reduce nullspace conditioning; the current full-field divergence residual
   is measured, not assumed zero. Independently evolve small perturbations and
   refine mesh/amplitude/time to distinguish physical growth from spatial error.
3. Fix only a demonstrated numerical defect. Do not suppress growing modes by
   dropping species feedback, fitting E/B/pressure, weakening residuals, adding
   artificial mass or repairing heat/energy. A stable uniform spectrum or fixed
   current species test does not certify the coupled heterogeneous state.
4. Repeat full contact and other unforced parent trajectories; close local
   enthalpy/Poynting interfaces. Then continue heated moving exact vacuum,
   transport/ionization, actual ends/coils, converged formation and experiment UI.
   Keep frc-7mm.4/.6/.7, frc-7mm and frc-goal open until their actual gates pass.

All previous source archives remain intact. No remote or publication; preserve
the user's modified outputs/frc_showcase/poster.png. Older sections below are
historical. This checkpoint made substantial new progress; it is not a blocked
goal and must not be resumed by resetting its scope or rerunning old remote jobs.

---

# Previous coupled moving-plasma checkpoint

The user made the full experimental-sandbox objective an active goal. Read the
new top section of GOAL_STATUS.md, REGULAR_PLASMA_NUMERICS.md and
REGULAR_PLASMA_VERIFICATION.md. The older isolated-species entry below is now
historical. The parent goal and eleven unfinished beads remain open.
**425 tests pass in 219.19 seconds.** Every job from this continuation is terminal;
there is no local or remote simulation to resume. Continue from saved evidence.

## Current evidence

- Actual fluid/current/two-species/magnetic RHS and full simultaneous RK4:
  25 frozen cases, complete bitwise replay, independent quadrature 16/20.
  Source 19f68fb08f949c4526bc1347d18652fb1eeb22d720d483a2074d6a71f7f40358.
- Seven full fast-wave periods; time-refinement ratios 15.999437/15.999863.
  Finest B/J signal errors 2.35089e-6/2.59599e-6; density/pressure 5.6764e-4.
  Contact has convergent original instantaneous rates only, not a full transit.
- Seven independent stable physical-energy audits, replayed at both quadratures;
  six slow-branch amplitude controls retain a weakest-pressure precision floor.
- The old full coordinate-matrix attempt's rejection was captured/replayed.
  **It is a localized near-axis toroidal magnetic perturbation, not the initial
  pressure-only hypothesis.** Inner physical and outer coefficient norms differ.
- Opt-in RefinedRegularPlasma tightens only the inner induction solve until all
  original unchanged outer residuals pass. Two full 360-coordinate matrices
  replay; independent linear errors 2.04814e-9/2.20696e-11. A complete coarse wave
  exactly matches the original. Source bc562813c4df25bb032e693b4e59adf0782aa9e92ee29f1233a36b34a1a2cd3f.
- Seventeen earlier calculation closures remain unchanged. All archives include
  complete source, arrays and metrics; do not overwrite or mutate frozen sources.

## Immediate next action

1. Read companions for regular_plasma, regular_plasma_step,
   refined_kernel_induction and refined_regular_plasma. Use the **opt-in refined
   class** for new work; original RegularPlasma remains the frozen baseline.
2. Complete a full exact FiniteMassContact transit on actual moving state with
   physical masses, both temperatures and full E. Derive a stiff time method or
   justify the explicit fast-mode cost; the contact transit is 100 microseconds,
   far longer than the 116.5 nanosecond fast-wave fixture. No rate clipping,
   artificial mass, field/pressure fit, heat or energy correction.
3. Preserve strict original residuals and atomic rejection across all stages;
   refine time/mesh and measure point/sampled positivity, initial/final physical
   norms and grid versus continuous conservation. frc-7mm.6 remains open until
   broader moving robustness/adoption is supported.
4. Close the original local enthalpy/Poynting budget with reconstructed face
   jumps. Integrated energy and nonlocal entropy graph work do not prove it.
   Inspect constant-potential gauge/EOS difference precision for very weak slow
   pressure signals; amplitude 0.001 controls retain a measured accuracy floor.
5. Only then continue heated moving **exact vacuum**, physical transport/
   ionization, actual coil/end boundaries, converged formation and experiment UI.
   Keep the original full goal open, including stability/confinement limitations.

There is no remote or publication. Preserve the user's modified
outputs/frc_showcase/poster.png. Existing study figures are verification reports,
not a proposed laboratory experiment. No old remote job requires continuation.

---

# Previous regular species checkpoint

The user resumed on 2026-10-04 and supplied a provisional lab envelope; research
is active and all ten previously unfinished beads remain open. **416 tests pass**.
Read the latest section of GOAL_STATUS.md, REGULAR_SPECIES_NUMERICS.md and
REGULAR_SPECIES_VERIFICATION.md. The species component is now implemented;
the earlier instruction below to begin that component is historical.

Frozen source `45a71e5b59cbcab6bb8556d531510d68da8389b16abcc82e1e61c0bc3237b78d`
has 24 bitwise-replayed cases, independent continuous audits, a local enthalpy
audit and eight full-transit fixed-current entropy trajectories. Actual moving
fluid/current evolution has not yet adopted it. Fifteen earlier calculation
sources are unchanged, and every job from this continuation is terminal.

## Immediate next action — frc-7mm.4 remains open

1. Read companions for regular_species, species_projection,
   regular_species_transport, regular_coordinate_fluid, regular_current,
   current_transport and kernel_induction before assembling a new coupled module.
2. Use the actual n,q,l,mz,Ce,Ci state and one common density derivative. Add
   species fluid_rate to bulk/current-stress rates; add **reduced_forcing** to
   induction, restoring **potential_electric** afterward for the full original
   electric/Ohm diagnostic. Do not numerically curl a large compatible gradient.
3. Verify independent original parent rates and complete unforced compressible
   wave/contact trajectories, retaining physical masses, local temperatures,
   axis/wall/max norms, time/mesh refinement and atomic failure. Derive the
   coupled time-step/implicit method for the actual fast electron/current modes;
   fixed-current RK4 stability does not select a safe coupled time step.
4. Close the local physical enthalpy/Poynting interface budget. The new graph
   transport is conservative and globally supported; it does not prove exact
   geometric face currents. The piecewise strong enthalpy divergence retains
   distributional jumps, and continuous energy remains separate from grid work.
5. Retain original fluid accuracy/spectral limitations and proceed to heated
   exact vacuum, physical closures, actual coil/end conditions and converged
   formation only after the coupled numerical gates pass. Keep frc-goal open.

The user's prospective 1 m quartz tube, 10–20 cm diameter, selectable gas,
0.1–100 mTorr neutral fill and adjustable Helmholtz-like bias are recorded in
RESEARCH_QUESTIONS.md. The verification fixture is not a selected lab experiment.
There is no Git remote. Preserve the user's modified outputs/frc_showcase/poster.png.

---

# Historical new-thread handoff

The user stopped this thread on 2026-10-04. The complete portable handoff is
[THREAD_HANDOFF.md](THREAD_HANDOFF.md), including all ten unfinished beads,
original scope, verified checkpoint, execution details and remaining gates.
No species implementation was started after research commit `7b3c8b2`.

The previous continuation note below is retained as historical context.

---

# Continue from the 2026-10-04 checkpoint

The user resumed on October 4. The overarching FRC goal remains open; this is
verified component progress, not a completed formation simulator. This repository
has no configured remote. Research and issue state are saved locally.

## Verified now

- Current advection, drift stress, stretching and self-drift: 16 frozen cases,
  bitwise replay, independent order-16/20 continuous-field audit. **408 full
  tests passed in 228.47 seconds.** Source `2f682eece854d54b3b0d33026709f8d34617cc4633d8b3875a0f3e2b8e0d9540`.
  The study CLI is complete; the missing records/audit modules from the prior
  pause have been implemented. frc-7mm.3 is closed. Read CURRENT_TRANSPORT_NUMERICS
  and CURRENT_TRANSPORT_VERIFICATION for exact scopes and measured limits.
- The old remote source-3f nr32 contact finished all 4116 steps and passed its
  original-platform final-step/physical audit. All 14 face-radial trajectories
  are complete and audited; frc-5b5 is closed. Evidence is in
  `outputs/face_radial_remote_completion/`. No restart occurred.
- Four species-pressure/entropy charge-compatibility probes, both quadrature
  orders replayed bitwise. Source c6a23e1a, `outputs/species_charge_probe/`.
  This diagnostic is separate from the 408-test solver baseline. It is not an
  implemented species transport operator.
- All 13 preceding calculation source closures are unchanged. Do not mutate
  their archived source definitions or rerun into existing immutable records.

## Next action: frc-7mm.4

Read `REGULAR_SPECIES_COUPLING_DESIGN.md` and the relevant module companions.
Derive the shared current/particle flux and both species entropy updates using
the actual regular n,q,l,mz,Ce,Ci state and physical h_e,h_i energy gradients.
The identity grad(p)/n=grad(gamma*h)-h*grad(s) exposes the required charge-flux
pairing. Compatible potential gradients annihilate current work, but their local
approximation error is nonzero and must converge; a global cancellation alone
is insufficient. Do not silently change density continuity or repair heat.

First implement and verify original local pressure and entropy rates with
nonuniform density, independently varying species temperatures and nonmatching
element traces. Measure axis/wall/max errors, species entropy, charge/particle
continuity, torque and physical versus grid energy. Then combine the verified
components and run unforced complete-parent wave/contact trajectories, including
enthalpy/Poynting boundaries and convergence. The lower-order continuous-energy
and earlier fluid order/spectral limitations remain recorded under frc-7mm.1.
frc-7mm.2 remains open for full coupling; frc-7mm.5 covers later adoption/formation.

## Later formation gates

Preserve exactly zero exterior plasma in heated moving-vacuum runs. Physical
transport/ionization, ends/coils and converged formation remain required.
Hardware stays abstract prescribed waveforms. Closed flux or axisymmetric
equilibrium does not establish 3-D stability or particle confinement.

## Machine and presentation

The authorized remote host is `m@192.168.0.202`; the completed old study is in
`/home/m/Code/frc-tool-comparison/face-radial-20260930`. The original PIDs were
absent on resumption and terminal evidence was audited/retrieved. No remaining
remote run is expected from that study; revalidate actual process identity
before any future action. Existing host authorization includes data transfer.
It was CPU NumPy/SciPy work, not a GPU-accelerated solver.

The six-panel convergence report is in
`outputs/current_transport_diagnostics/convergence.png` and `.pdf`.
The earlier interactive field presentation is preserved in
`outputs/frc_showcase_checkpoint/` (HTML, poster, MP4). It still depicts the
saved coupling case 07, with camera orbit only, and its historical status label
is from the previous checkpoint. Rebuild/serve via
`.venv/bin/python -m examples.showcase_checkpoint --serve --port 8779`.
The user's modified `outputs/frc_showcase/poster.png` remains outside our commits.
