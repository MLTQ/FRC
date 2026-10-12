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
