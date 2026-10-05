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
