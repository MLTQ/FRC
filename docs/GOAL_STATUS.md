# Paused objective: a credible experimental FRC simulation

## User-requested pause and handoff (2026-09-30)

Paused for next week to conserve credits. See [NEXT_SESSION.md](NEXT_SESSION.md)
for the exact restart point, ordered beads and remote-study observation.
The overarching goal remains open. No new research study was launched for
this handoff; no remote job was stopped or restarted.

Beads **frc-7mm.3**, **frc-7mm.4**, **frc-7mm.5** track current-transport audit,
remaining species coupling and moving-vacuum/formation validation respectively.
The saved transport candidate passes **4 focused tests in 1.66s**; it has no
full regression or frozen 16-case audit yet. Its study CLI still requires the
records/provenance module. Keep it distinct from the **402-test / 16-case**
verified coupling baseline below. No production solver was replaced.

The updated viewer, poster and camera-orbit movie are in
`outputs/frc_showcase_checkpoint/`, using frozen verified coupling case 07.
They show actual saved vector fields in a filled periodic cylinder, not a
new formation result. Camera motion is not physical time evolution.


## Last fully verified checkpoint: physical regular-fluid/current bridge (2026-09-30)

The previous goal turn was verified progress: it preserved the failed Gauss-r
representation and verified regular-axis mixed fluid transport. This turn
connects that actual point density and physical velocity to compatible current
and magnetic work. No cell-mean reinterpretation or alternate continuity was
introduced. Force work, actual integrated axial force and angular torque use
the derived physical adjoints. Physical current inertia and its density-energy
pullback now use the same actual regular density and continuity derivative.

The full suite passes 402 tests in183.03s. Sixteen degree/resolution/quadrature/
density cases pass original-platform bitwise array/metric replay and independent
root-polynomial/FFT physical integration. Frozen source is
8c9339a6eaf0c88b88bffe108d5ca2c535a8bfbf33fa4bd11b0c0aaebba6c784.
All twelve previous source closures are unchanged. All new local study, test,
audit and reporting jobs are terminal. The six-panel report was inspected.

Degree5 smooth-force errors at2/4/8/16 fluid elements decrease to radial
9.10329e-10, azimuthal4.26004e-10 and axial8.67456e-12; first-element errors
reach roundoff/quadrature limits. Force work, integrated force and torque agree
within5.75e-16. Actual Ampere/Lorentz force is separately checked: finest global
component errors7.44700e-8,2.06750e-7,1.04214e-6; first-element axial1.38035e-6.
Manufactured magnetic-rate errors at8/16/32 field cells decrease4.62633e-4,
1.58115e-5,5.05377e-7. Separate field/fluid-grid controls expose both sources
of error. Actual continuity replaces analytic n_t in changing-density inertia.

The ideal force/electric pair agrees within9.18e-16; Hall work within7.59e-17.
The current-energy density pullback matches complex differentiation within
4.36e-16. Original induction residuals remain below1e-10. Independent physical
energy/work integration agrees closely; its largest relative discrepancy is
3.68e-13 for a nearly cancelling coupled magnetic integral. See the full audit
and REGULAR_COUPLING_NUMERICS/VERIFICATION for exact scopes and norms.

IMPORTANT: this remains an instantaneous, incomplete moving pair. Its total
energy rate contains integral mu*J²*n_t/(2e²*n²). At the finest matched fixture,
measured total rate is−0.001498461177W versus derived exchange−0.001498461105W.
The missing current-advection/drift-stress transport is explicitly retained,
not relabeled heat or canceled with a force/energy repair. No complete coupled
time evolution or formation claim follows from this component's successful audit.

The older frozen remote3f study remains13/14 complete and audited. At
2026-09-30T18:06:01.201528Z the nr32 contact is live at2050/4116 accepted steps.
Parent147751 and workers147839/147845/147852 are still in the correct
/home/m/Code/frc-tool-comparison/face-radial-20260930 directory on
m@192.168.0.202. Saved status is regular_coupling_diagnostics/remote_status.json.
Never restart from an observation timeout; no GPU acceleration is claimed.

NEXT continue frc-7mm.2: derive the shared current-advection and drift-stress
work using the actual regular continuity and cylindrical vector derivatives.
The continuum cancellation identities are now explicit in REGULAR_COUPLING_NUMERICS
and a candidate weak advection identity in REGULAR_CURRENT_TRANSPORT_DESIGN;
a global fitted cancellation is insufficient. Complete physical species
pressure/entropy and enthalpy/Poynting transport, then original parent-rate
and unforced wave/contact evolution controls before adoption. Prior unforced
fluid order and degree3 spectral limits remain recorded under frc-7mm.1.
Heated moving exact vacuum, physical closure/ionization, real ends/coils,
converged formation and3-D stability/confinement remain open. No production
solver was replaced. Full goal remains open and is now paused; state is local without a remote
Git repository or publication.

## Previous checkpoint: regular-axis mixed fluid transport (2026-09-30)

This continuation found and resolved a physical representation defect in the
previous fluid prototype. Mixed radial/axial/rotating flow revealed large
between-point velocity spikes despite conserved grid energy. At degree 5,
coarse old unforced grid/continuous energy gaps reached 0.06675, while a
physical axis-metric eigenvalue reached 576461. A naive physical force bridge
also failed the analytic f(r)=r test with nonconverging relative error 0.66543.
The old solver and all fourteen comparison trajectories remain preserved.

The new regular-coordinate state uses x=r²/2 and actual n,q=n*vr/r,
l=n*vtheta/r,mz,Ce,Ci. True angular momentum is 2*x*l. Axis velocities are
dependent regular traces with no zero-energy independent axis velocity slot.
Derived geometric work and joint entropy-density fluxes conserve grid physical
energy and angular momentum. The constant-density first-element velocity
metric is also the exact continuous physical metric. No energy reset, force
fit, artificial viscosity, radius/density floor or atmosphere was introduced.

The full suite passes 398 tests in 181.18 s. Eight original-equation rate
cases, eight complete acoustic spectra and fourteen mixed trajectories have
source/raw hashes and bitwise/independent physical audits: thirty new records,
plus fourteen preserved old comparison trajectories and three metric probes.
New source is 06dfff4762b721a154d9b8522f935427da0a13f335ac74a69549a08cdc958e7e;
old mixed source is a3d0b0c387bd4696fca958c20038d41ad7725746f4809d4052f80c36493e07b6.
Ten earlier source closures and both current closures are verified unchanged.

At degree 5 and 2/4/8/16 radial elements, manufactured continuous density
errors decrease 1.82947e-5,3.62103e-7,7.05634e-9,1.64671e-10. Fixed central
velocity RMS/c errors decrease 9.69768e-6 to 4.29509e-11. Unforced sampled
speeds remain near 0.29646*c on every mesh; final grid/continuous energy gaps
refine 5.72337e-8,5.16148e-9,3.63098e-10,1.70248e-11. Independent higher
physical quadrature agrees within 2.76914e-16 across all saved new snapshots.

Important limits remain measured, not concealed: unforced successive-mesh
velocity differences converge only roughly second order, and the degree-3
third sorted positive acoustic mode is not the third Bessel root. Higher
manufactured order does not resolve either issue. Grid conservation is not
exact continuous variable-density energy/angular conservation. The prior
part of the accumulated manufactured-work ledger is not independently
replayed; the last increment and every saved physical history are audited.
See REGULAR_COORDINATE_NUMERICS/VERIFICATION and the inspected six-panel
report in outputs/regular_euler_diagnostics. All new local jobs are terminal.

The older frozen remote3f study remains 13/14 complete and audited. At
2026-09-30T17:36:47.990170Z, its nr32 contact checkpoint is 1435/4116 accepted
steps. Parent147751 and workers147839/147845/147852 remain in the original
/home/m/Code/frc-tool-comparison/face-radial-20260930 directory on
m@192.168.0.202. Status is saved in regular_euler_diagnostics/remote_status.json.
Do not restart from an observation timeout. No GPU acceleration is claimed.

NEXT frc-7mm.2: derive and independently test the current/magnetic bridge on
the actual regular-coordinate density and velocity fields, with physical
variable-current inertia, changing-density work and true angular torque.
Start with arbitrary smooth force/torque and variable-density pullbacks; do
not silently map the new point state to the old cell/face representation.
frc-7mm.1 stays open for remaining transport limits; frc-5b5 awaits the final
remote audit. Full stress and enthalpy/Poynting transport, heated moving
exact-vacuum edge, physical closure/ionization, real boundaries, converged
formation and 3-D stability/confinement remain unproved. No coupled production
solver has been replaced. This is verified progress toward an active goal.
All repository and issue state stays local; there is no Git remote/publication.

## Previous checkpoint: cylindrical energy and angular transport (2026-09-30)

This continuation makes verified component progress; the full FRC goal remains
active. A separate cylindrical Gauss/SBP point state now jointly transports
density, physical radial momentum, **actual angular momentum**, axial momentum
and both species entropy densities. Its derived centrifugal/pressure work pairs
with actual physical kinetic/thermal energy in the declared grid quadrature.
No energy reset, force fitting, viscosity, density/radius floor or atmosphere
was introduced. Existing coupled plasma solvers were not silently replaced.

The full suite passes391 tests in165.29s; six focused tests pass again after a
reporting-only correction. Seventeen full2-D original Euler rate cases, six
independent Bessel/acoustic spectra and seventeen unforced trajectories pass
source/raw hashes and original-platform bitwise/independent physical audits.
New frozen source is14223ded4171ee70b3b790323913b1ff302d994df98a13c6328522273ebcf6ef.
All nine earlier source closures are verified unchanged.

Degree5 first-element radial acceleration errors decrease3.09383e-4,8.64907e-6,
1.24677e-6,9.05808e-8,5.87185e-9 at2/4/8/16/32 radial elements. Finest maximum
pointwise relative axis error remains2.63e-7, separately exposed. Full errors
reach9.92e-10. Degree7's fixed128-axial floor is confirmed by256-point controls,
which reach1.29e-11; its finest axis errors become nonmonotonic near1e-11.
These element counts are not equivalent to earlier face-momentum mesh counts.

Grid physical-energy rate defects are <=6.68e-15; state-rate conservation
<=8.74e-15. Each radial element's energy balances its interface flux in tests.
Full-contact density errors decrease2.37e-4/4.26e-6/6.92e-8 at16/32/64 axial
points. Fixed32-point energy drift decreases−1.50030e-11/−9.30478e-13/−5.79737e-14
under step halving. History invariants conserve within3.39e-15. Actual nonlinear
radial waves complete a full period; amplitude controls distinguish physical
nonlinear departure from their linearized Bessel target. Rotation runs span one
radial sound-crossing, not a complete rotation or plasma stability experiment.

Independent continuous fields use even polynomials in r² for n,mz,Ce,Ci,mr/r,
L/r², so the axis is regular and kinetic energy finite. This differs from the
SBP differentiation polynomial in r. Grid conservation is not exact continuous
energy/angular conservation. Across saved trajectory snapshots, higher physical
quadrature agrees within4.35e-16, while grid/continuous energy differs by up to
8.83e-11. Rotating continuous angular drift refines1.40e-9→3.71e-12→1.84e-13.
All measured differences remain in the evidence, not repaired away.

The first reporting attempt is preserved in cylindrical_entropy_before_velocity_norm:
relative velocity error was undefined at a zero-velocity turning point. Correct
RMS/sound-speed normalization changes no solver array;153 shared raw files replay
bitwise. Final evidence and visually inspected six-panel figure are in
outputs/cylindrical_entropy_verification and cylindrical_entropy_diagnostics.
See CYLINDRICAL_ENTROPY_NUMERICS/VERIFICATION for equations and limitations.

The older remote3f nr16 contact has completed2060 steps, passed its original-
platform audit and been retrieved. Together with prior local/remote cases,
13 of14 are complete/audited. Its continuous density/Te errors are1.918e-7/
8.584e-7. Only nr32 contact remains live: at2026-09-30T16:47:37Z its accepted
checkpoint is820/4116, parentPID147751 and workers147839/147845/147852 remain
in /home/m/Code/frc-tool-comparison/face-radial-20260930 on m@192.168.0.202.
See outputs/face_radial_remote_controls/remote_status.json. Do not restart a
live worker from an observation timeout; no GPU acceleration is claimed.

NEXT frc-7mm.1: verify fully mixed nonlinear r/z evolution with mesh/time
controls, then derive the variable-current/magnetic bridge on this actual
common state. No unproved nodal-to-cell/face conversion may be adopted.
Full species/Maxwell/drift stresses and enthalpy/Poynting work, heated free-edge
exact vacuum, physical transport/ionization, real ends/coils, shocks/kinetic
validity and converged formation all remain open. Audit the one remaining
remote trajectory when terminal without changing its source. Work/issue state
is local, with no Git remote or publication. This is verified progress, not a
completed goal or external blocker.

## Previous checkpoint: joint physical energy transport (2026-09-30)

The previous goal turn was verified progress: common radial momentum corrected
the near-axis force mismatch and preserved the incomplete energy contract.
This continuation derives a joint mass/momentum/two-species entropy-density
flux that conserves actual physical energy in its declared nodal integration
norm. It adds no energy reset or force-only correction. This is a new point
state, not a silent reinterpretation of the existing cell/face state.

The complete suite passes **385 tests in 164.90 s**. All 12 original axial Euler
rate cases and 13 unforced full-transit/wave trajectories are independently
audited, with bitwise state/rate/last-step/history replay. Independent 8x/16x
physical integration agrees within 1.92e-16; nodal-versus-continuous physical
energy discrepancy is separately retained, at most 6.70e-14 in saved states.
The four-panel convergence/energy figure was visually inspected.

At 128 points and order6, physical-energy-rate relative error is 6.01370e-9;
electron/ion pressure-rate errors are 4.01490e-10/1.19765e-8. Independent local
H_t matches conservative energy-flux divergence within 2.13e-14; its global
rate defect is at most 1.02e-16. Actual nodal mass, three Cartesian momenta and
both entropy-density histories conserve within 1.79e-15. This is not yet a
cylindrical angular-momentum theorem.

Finest contact/wave continuous density errors are 1.09227e-9/3.15956e-10;
their density-perturbation errors are 5.26361e-9/4.48378e-9. At fixed 32 points,
contact energy drift decreases from −2.46222e-11 to −1.52514e-12 to −9.40984e-14
under step halving, approximately fourth order. Continuous physical-energy
integration confirms the same trend. Reverse/time controls, separate species
errors and initial representation differences remain in the records.
Global entropy conservation does not give exact local isentropy: the analytic
isentropic wave acquires specific-entropy errors that refine from 1.15e-6 at
16 points to about 9.96e-12 at 128. This limitation is explicitly retained.

New frozen source is
924a586a84df69c28fe1493bf73e285efcecc02e63b6825926f5be232f94170d in
outputs/two_entropy_verification. Equations, state contracts and remaining
scope are in TWO_ENTROPY_NUMERICS/VERIFICATION. Existing source closures and
the running 3f4120fa radial prototype are unchanged. All new local computation,
test, audit and report jobs are terminal. No production plasma adoption occurred.

Seven additional remote radial-face controls completed and passed their
original-platform audits, then were retrieved to outputs/face_radial_remote_controls.
Combined with the five local cases, 12 of 14 are complete/audited. Finest nr32
acoustic n/Te errors are 1.69024e-8/1.12770e-8. The nr8 quarter-step contact
retains −4.10007e-9 physical energy drift with only −2.53e-15 finite-step
remainder, confirming that separate prototype's spatial transport defect.

Remote worker state is verified at 2026-09-30T16:00:33Z: parent PID147751 and
workers147839/147845/147852 remain in the correct study directory
/home/m/Code/frc-tool-comparison/face-radial-20260930. Saved accepted checkpoints
are nr16 contact721/2060 and nr32 contact205/4116. The other assigned cases
are complete. See outputs/face_radial_remote_controls/remote_status.json.
Never restart from an observation timeout; use existing checkpoints only after
authoritative evidence that a worker stopped. No GPU speedup is claimed.

NEXT frc-7mm.1 must extend joint transport to the actual cylindrical state:
true radial kinetic/thermal metric, regular-axis force, conservative angular
momentum, geometric torque and impermeable-wall response. Verify the same full
two-dimensional original equations and radial modes before adoption. Then
derive/verify the common variable-density current/magnetic and full species
enthalpy/Poynting work; no unproved nodal-to-cell/face bridge is allowed. Audit
the two remaining remote contacts when terminal without editing their source.
Full heated-edge/exact-vacuum, physical closure/ionization, real boundaries,
converged FRC formation and experiment requirements remain open. This is
verified component progress, not completion or an external blocker.

## Previous checkpoint: common radial physical momentum (2026-09-30)

This continuation made verified component progress; the full FRC goal remains
active. The complete regression suite passes **382 tests in 169.23 s**. Nine
original two-dimensional Euler rate cases and five completed unforced local
trajectories pass source/raw-array, bitwise replay and independent physical
energy/history audits. Seven radial acoustic spectra and eight independent
polynomial metric probes are retained. The six-panel scientific figure was
visually inspected, including corrected mesh-axis labels.

Independent interior radial faces now define the actual momentum field used by
continuity, velocity and the physical kinetic metric. The old mixed cell/flux
adjoint was inconsistent; replacing its metric alone would have been singular.
The common radial-face method reduces degree5 first-annulus acceleration errors
from the old4f sequence4.956/6.330/6.651/6.732 to
1.75954e-5/1.91939e-7/3.08493e-9/1.19443e-10 atnr8/16/32/64.
Whole-domain acceleration error reaches1.32683e-9. Actual pressure, centrifugal
and advective components are separately checked near the axis. Their individual
outer-wall constrained responses are not mislabeled as unrestricted forces.

The intermediate exact radial total-work correction is preserved with its
failed convergence test. The current prototype evolves physical centrifugal
force and exact pressure/thermal work pairing, while retaining the convective
energy discrepancy. Its scaled global energy-rate defect decreases from
8.28819e-7 to3.51141e-12 under refinement but is not exactly zero. The nr8
full-transit contact completes with n/Te errors2.90525e-5/5.15852e-5, physical
energy drift−4.10015e-9, account drift4.35e-16 and finite-step remainder−7.41e-14.
Acoustic and rotation nr8/16 also complete. No energy or temperature repair,
artificial viscosity, density floor or exterior atmosphere is introduced.

Source3f4120faf1b96ee2e94453339efe215399071fd632827aa6ae459cd9c5126737 is
frozen in outputs/face_radial_rates and outputs/face_radial_pilot. The new
equations, limits and measured results are in FACE_RADIAL_NUMERICS and
FACE_RADIAL_VERIFICATION. Existing verified source closures remain unchanged.
The prototype has not replaced the coupled moving plasma solver.

Remote state checked at2026-09-30T15:34:58Z: detached parentPID147751 and
workers147839/147845/147852 are running in
/home/m/Code/frc-tool-comparison/face-radial-20260930 on m@192.168.0.202.
Nine remaining protocol cases5,8,9,10,11,12,13,1,2 are assigned there. At the
check, nr32 rotation was complete but unaudited; accepted checkpoints existed
for nr32 acoustic90/208, nr8 quarter-step contact618/2060 and nr8 reverse
contact510/1032. See outputs/face_radial_diagnostics/remote_status.json for
the observation, not a claim that all requested cases have finished.
The remote protocol uses its CPU/SciPy backend; no GPU acceleration is claimed.

Long runs now save actual accepted state/history/energy checkpoints about every
5percent. A nonlinear interrupted/restarted wave reproduces uninterrupted arrays
bitwise. Resume only after proving the worker stopped; use --resume on its
existing directory and original platform/source. Do not restart a live job or
rerun terminal failed cases merely because an observation times out.
All local numerical/test/audit/report jobs are terminal at this checkpoint.

The previous remote4f worker59269 is terminal (SSH exit255/reset); subsequent
/proc inspection found no Python worker in its original directory. Thirteen
saved cases were audited on the original host and retrieved to
outputs/radial_work_remote_recovered:12complete,1failed with its actual rejected
step replayed. Case2, nr32contact, is explicitly interrupted/missing after a
last log at75percent and7198.1s; it is not counted complete or silently rerun.
No further hours are spent on this superseded force-inaccurate formulation.

NEXT finish/audit the remaining checkpointed controls without editing their
source. frc-5b5 remains open pending those trajectories. New frc-7mm.1 tracks
joint conservative physical convective kinetic-energy/momentum transport;
derive it without the already-disproved force-only translation correction.
Then integrate full magnetic/current/Maxwell/drift and enthalpy/Poynting work.
Heated free-edge exact vacuum, physical closures/ionization, real boundaries,
converged formation and the full FRC objective remain open. Repository and
issue state are local; no Git remote or public publication is introduced.

## Previous self-consistent fluid checkpoint (2026-09-29)

The new smooth fluid components now advance actual density, all momenta/angular
momentum and both physical log entropies. True reconstructed kinetic energy and
its changing-density work remain separate from a conservative energy account.
The latest complete regression run passes **376 tests in 151.46 s**.
No old verified current/thermal/bulk source closure was modified.

Longer tests exposed a radial acoustic instability in the first conservative
pressure stencil. Baseline source 0aa417e5 is preserved, including failed contact
states. Seven spectra with original-RHS difference and matrix-exponential/time
controls prove spatial growth. Nine uniform-translation controls also prove
that a pressure-force-only exact-energy correction would add net momentum.

Separate radial thermal/kinetic work pairing removes the growing linear modes
and strongly improves rotating force balance. Its nine 2-D parent-rate cases
and independent component audit nevertheless reveal a nonconverging relative
radial acceleration error in the first annulus (about 6.7 at nr64). Thus it is
experimental and not adopted in the coupled plasma solver. See LOG_FLUID_NUMERICS
and RADIAL_WORK_NUMERICS. Candidate source 4f0c824a662e581dc78ddfbd81a2c9cc06342fb85d19d62fda44d0fb709ab3a9.

The revised nr8 contact completes a full transit at half the original step:
density/Te relative L2 2.90525e-5/5.15852e-5, physical energy drift −4.10015e-9,
conservative account drift −1.45e-16 and separated finite-step energy defect
−7.36e-14. The quarter-step run also completes: final density/Te differences from
half step are6.29e-11/8.79e-11; remaining finite-step energy defect is−2.67e-15.
The original larger-step failure remains archived. Eight local trajectory records
and nine parent rates are independently audited. The inspected six-panel figure
and LOG_FLUID_VERIFICATION report retain the local radial-force limitation.

At that checkpoint remote candidate59269 was running at
/home/m/Code/frc-tool-comparison/log-fluid-20260929. It has since stopped;
the current recovery record above supersedes that live status. Baseline15653 is terminal
(exit1: five failed contacts, nine complete cases). All14 records are audited on
their original platform and retrieved to outputs/log_fluid_remote_baseline.
All five actual rejected steps replay with unchanged accepted arrays. Their
already unstable last states expose quadrature errors up to1.05e-3, with higher
quadrature still changing the worst by1.73e-4; this warning and the original audit
rejection are preserved, not marked converged. Both audit workers50776/40535 are
terminal. All local numerical/audit/report jobs are confirmed terminal.
Original preview11395 remains available at port8778, separate from computation.
No remote exists for Git. User-approved SSH data transfer has resumed; no public
publication occurred. Mac disk pressure was relieved by removing only66.97MB
of regenerable .venv bytecode, with about410MB free at the last check.

Four additional force-term probes isolate nr64 first-annulus pressure and
kinetic-source relative errors6.44600/6.44732, while convection is2.15e-10.
Track correction as frc-5b5; balanced rotating forces can hide these compensating
errors. The five previous verification source identities remain unchanged.

The next action at that checkpoint was to finish/audit the fixed-source
comparison and resolve the near-axis radial pressure/kinetic work adjoints before
full pressure/Maxwell/drift and current/enthalpy/Poynting coupling. The heated
free edge/exact vacuum, physical closures, formation and full FRC goal stay open.

The goal is a numerically verified and physically justified simulation in which
we can vary geometry, plasma and prescribed coil waveforms to investigate FRC
formation. Passing the existing stationary-conductor tests is not completion.

Latest comparison checkpoint: **17 WarpX GPU cases audited**, with a working
isolated RTX 4090 installation and a source-level Gkeyll fit review. Retain our
existing bench and custom axisymmetric path; an external engine is not yet a
verified shortcut around the physical plasma-vacuum edge. Full-Maxwell vacuum
waves converge at order two. The tested hybrid waves retain measured time/energy
errors; no FRC formation or GPU speedup claim is made. See ENGINE_COMPARISON.
No production backend or Stage 2 gate changed.

Latest custom checkpoint: **370 passing tests (126.72 s), 14 bulk thermal
rates and 10 changing-density trajectories independently audited**, plus six
bitwise-replayed diagnostic pressure probes and 24 linear face controls.
One conservative particle flux now advances density and both physical species
log entropies under prescribed bulk flow. The initial unstable normal-face fit
was diagnosed, preserved and corrected with centered faces, without added
viscosity or density floors. Finest full-transit density/electron-temperature
errors are 8.34e-9/6.57e-9; finest compression errors are 1.37e-6/8.27e-7.
Physical work and entropy identities, finite-step drift and initial projection
errors remain distinct. See BULK_THERMAL_NUMERICS/VERIFICATION and the inspected
six-panel report. Source a74e877fe9dd2030a6957d717d87dc1a1e96c4866fdb92ba29e5b76c157fcd2b.
All scientific workers and audits are terminal. The storage-error interruption
was recovered by recomputing only five missing cases under unchanged source.

The user explicitly authorized transfer to the GPU host after automatic review
initially rejected it. A first copy of 956,639,669 bytes of source, tests, docs,
issue export, configurations and outputs completed successfully at
m@192.168.0.202:/home/m/Code/frc-tool-comparison/recovery-bulk-a74e877-20260929/.
No Git remote was configured and no research was publicly published.
The user-requested 3D showcase is complete: an offline interactive viewer and
a visually inspected 2800×1900 current/flux-density poster from source_g1_32.
All three saved times, fixed scales, synchronized rotation and export were
checked. See SHOWCASE_3D. This initialized resistive FRC fixture is separate
from the new bulk thermal component. Only its localhost preview remains live
(session 11395, port 8778); it is not a scientific worker.

Full nonlinear adoption is still incomplete. Next derive conservative momentum
and angular-momentum evolution with the true physical kinetic metric,
velocity/vorticity, pressure/Maxwell/drift stresses and complete
current/enthalpy/Poynting work, then independent parent rates and full nonlinear
wave/contact trajectories. The diagnostic pressure force converges but is not
an exact thermal-work adjoint; it has not been adopted. The earlier smooth
solver retains fast-wave pi error 5.38315%, manufactured nr32 pi-change error
70.2663%, and slow nr64 B_t error 16.8431. Previous log-thermal, variable-current,
compatible-wave and smooth-moving source identities remain unchanged.
Heated-edge/exact-vacuum, physical transport, formation and accurate interactive
experiment requirements remain open. This is progress, not goal completion.

## Completion evidence still required

1. Conservative axisymmetric compressible dynamics: mass, momentum, pressure/
   total energy, magnetic advection and resistivity, with compatible div B.
2. Verified plasma–vacuum coupling and suitable chamber/end/outer boundaries;
   no undisclosed conducting or massive numerical atmosphere.
3. A defensible transport/ionization model for the selected regime, with local
   scale/collisionality checks; Hall/two-fluid or kinetic treatment if required.
4. Independent tests of waves, shocks, magnetic forces, diffusion/heating and
   equilibrium preservation, including convergence and energy budgets.
5. Converged formation experiments with separate reversal, topology, force
   balance, lifetime and model-validity criteria. Robustness to time step, mesh,
   boundary and transport choices; no claim of stability from axisymmetry alone.
6. Scriptable and interactive experiments with complete reproducibility records
   and diagnostics that do not overstate the supported conclusions.

## Current evidence

The initial turn was progress: Stage 0/1 were implemented with 33 passing tests
and three waveforms plus a 27-case sweep. Their topology durations are not
converged and their illustrative argon parameters raise missing-physics warnings.

The first continuation implemented a conservative dynamic r-z MHD kernel that
passes independent cylindrical fast-wave, 2-D flux-advection, exact Sod shock,
resistive-decay/heating, axis-balance, divergence and temporal-order checks.
The next continuation added limited second-order reconstruction of fluid and
staggered magnetic traces, exact regular-axis moments, total-pressure source
integration and improved one-sided wall gradients. Smooth wave/advection errors
are much smaller; axis force accuracy improves, and a particular straight-field
polynomial equilibrium is held to roundoff. First-order comparison remains
selectable and is rerun on identical grids.
The bounded-reference checkpoint had **62 passing tests**. At identical resolution the
fast-wave error is 216 times smaller and the 2-D flux-advection error 129 times
smaller than the first-order reference on the reported benchmark grids.

An independent continuum Grad–Shafranov fixture now tests genuinely curved,
closed poloidal flux with reversed Bz. It starts in equilibrium, not from coil
formation. Its artificial motion decreases with refinement, but it is not held
exactly. A lower-pressure p0=0.1 Pa stress case previously failed reconstruction
on coarse meshes despite positive cell pressure. That failure is now fixed by
local admissibility reconstruction, with consistent geometric-source integration.
Cell averages and magnetic traces are unchanged by this limiter; there are no
cell floors or energy repairs. The original failure remains reproducible with
the explicit disable switch. Limiting affects 12.10% and 2.19% of accepted
cell-stages at nr=12/24, and zero at nr=48/96 in that fixture.
The default p0=1 Pa case completes across refinement after correcting the wall
derivative. Smooth temporal refinement verifies RK2; the quarter-period wave's
limiter-switching time ratios are separately retained and less regular.
For p0=1 Pa, the finest 192×384 equilibrium run has maximum artificial speed
0.00095046 times the reference Alfvén speed after 0.1 s, with flux relative L2
drift 2.1181e-5. The initial force volume RMS converges at approximately order
1.5, so this remains a measured limitation rather than exact equilibrium balance.
Fixed-grid time refinement barely changes that drift, confirming that spatial
error dominates this particular case. Mass/energy and face divergence remain
at floating-point accuracy. Generated figures were visually inspected.

Additional equilibrium stress cases down to p0=1e-5 Pa finish conservatively,
but their coarse-grid artificial motion is not physical accuracy evidence.
An independently solved magnetized double rarefaction remains positive and
conservative with unchanged uniform B. Its bulk L1 errors converge, while the
central pressure is still 114% too high at 640 axial cells. The admissibility
limiter is inactive there: this exposes underlying flux-discretization error,
not a pressure floor. It is tracked as an accuracy requirement, not hidden by
the passing positivity checks. No full-update positivity theorem is claimed.

Measured evidence is in `docs/MHD_VERIFICATION.md`,
`docs/EQUILIBRIUM_VERIFICATION.md`, `docs/POSITIVITY_VERIFICATION.md`, and their
JSON/plot outputs. The limiter derivation is `docs/POSITIVITY.md`. These remain
bounded verification solvers, not the coupled coil-driven formation experiment.

## Dynamic vacuum coupling checkpoint

That coupling checkpoint had **74 passing tests** (18.08 s).
`VacuumCoupledMHD` advances compressible gas inside a fixed, gas-filled chamber
with reflecting slip walls/ends and magnetically transparent boundaries. A
larger magnetic-only exterior satisfies the elliptic induced-current constraint;
there are no exterior density, pressure or conductivity degrees of freedom.
An adjoint magnetic force/advection pair, implicit midpoint and gas Euler Strang
splitting account for induced self-field energy, gas energy, prescribed-source
work and positive physical Ohmic heat. This response-energy budget excludes
external-field self/mutual energy and is not a coil-circuit energy estimate.
Failure leaves accepted state, time and cumulative work intact. A source time
cap is explicit because fluid CFL cannot discover arbitrary short pulses.

Independent compact-field J cross B errors decrease from 0.2707 to 0.006437 in
volume L2 as radial resolution increases 8→64, including axis cells. Passive
advection errors decrease 0.008325→0.002155→0.0005438 at 12/24/48 cells; Gaussian
resistive-diffusion errors decrease 0.007174→0.001822→0.0004573. Fixed-grid time
halving reduces momentum error by 4.010 and 4.050. A separate infinite-space
circular-loop integral distinguishes mesh error from magnetic-box error:
at 32 cells/m, zero-flux outer extents 2/4/8 m give 5.429%, 0.6791%, 0.08840%
probe error. With independent exact boundary data the 32-cell/m error is 0.03344%.

Two prescribed coils now induce current, gas acceleration and heat in this
reference. Across the reported dynamic cases the maximum relative response-
energy residual is 1.63e-14; mass change is at most 1.42e-16. Compatible face
divergence is near roundoff. **Pulse flow is not yet converged:** final peak speed
increases 0.0661→0.1473→0.1798→0.1999 m/s at 8/16/24/32 radial cells. The magnetic
boundary also matters: at 16 cells, extent 2→4 m changes final source work by
about 19%, and 3→4 m still changes it by about 3.8%. A small solver residual is
therefore not physical boundary accuracy. Figures were inspected; full inputs,
histories, source fingerprint and versions are retained with the report.

See `docs/COUPLED_NUMERICS.md`, `docs/COUPLED_VERIFICATION.md` and
`outputs/coupled_verification/metrics.json`. Earlier reports remain historical
measurements of their recorded source, not reruns with this different solver.

## Current development step

The separate fixed-material `MaterialCurrent` component now includes physical
electron inertia, spatially prescribed positive carrier density/nonnegative
resistivity, true magnetic vacuum and explicit electron drift kinetic energy.
Conditional Coulomb transport and local Hall/inertia/drift diagnostics are
integrated into reference experiment metadata. Prescribed conductivity and
moving dynamics are unchanged. The full suite has **192 passing tests (79.10 s)**.
All 35 current-component cases, an independent time audit, 43 conditional
states and 41 prescribed-edge samples are complete. Both final figures were
visually inspected; complete source/configuration/field hashes are checked.
No numerical job remains live at this checkpoint.

CURRENT_VERIFICATION.md separates continuum error, time order, physical current
energy and numerical false reversal. The unresolved de/R=.05, nr=8 screening
case produces negative Bz down to -0.0033032 of applied bias, despite a small
energy residual; the independent exact field is positive. Finer sampled fields
remove this artifact but nr=64 still has 7.3078% magnetic-field L2 error. The
inertial decay's nr=64 current/heat errors are 0.082887%/0.10486%; initial weak
surface-current error limits convergence and is retained, not hidden.

The original Ar+ reference has di/R=4.0972 and conditional electron mean-free-
path/R=0.11088. It does not justify ignoring Hall physics. Smaller physical
edge-gradient lengths raise electron inertia, carrier-drift and kinetic
requirements. The heated sharp-edge incompatibility in SURFACE_COMPATIBILITY
remains unresolved; the new linear fixed-density component does not fix it.

The next bounded coding step is the axisymmetric Hall/toroidal reference and
independent wave-dispersion verification in EXTENDED_PLASMA_PLAN.md. Its parent
species equations, current-energy force/transport terms, separate electron/ion
thermal evolution and justified edge treatment must precede complete moving
coupling. Spatial coil force/work checks remain in MOVING_COIL_PLAN.md. Actual
gas ends, mesh/shock robustness, converged formation and the experimental UI
remain required. The original low-pressure scope is retained; a denser/easier
fixture cannot silently replace it. `Experiment(stage=2)` stays gated and the
full objective remains active.

## Infinite magnetic exterior checkpoint

That checkpoint had **83 passing tests** (20.90 s). The preceding turn
made verified progress, with no external job remaining live at this checkpoint.
`vacuum_boundary="multipole"` now represents the induced magnetic field to
infinity using regular-axis decaying Legendre modes. Boundary nodal traces,
interior magnetic energy and an independently integrated exterior energy form
one symmetric matrix. Exterior nodes and modal amplitudes satisfy stationary
vacuum rows; modes have no fluid mass, current density or resistivity. The
midpoint response-energy budget includes energy beyond the field grid.
The historical fixed zero-flux boundary remains an explicit comparison.

Independent displaced toroidal-current quadrature excites both odd and even
modes. Interior flux relative L2 errors at 8/16/32/64 fluid radial cells are
0.01257, 0.0008744, 0.0001845 and 0.00004552. Interior/exterior normal-derivative
mismatch decreases 0.01653→0.003882→0.0009413→0.0002318, including matching-surface
corners. Exterior B error at fixed 32 cells drops from 32.3% with one mode to
about 0.0364% with twelve; higher modes then reach the spatial-error plateau.
A separate analytic dipole and independent all-space volume integral check
normalization, axis fields and the sign/magnitude of exterior magnetic energy.

All fixed/open pulse comparisons were rerun from current code. At 16 fluid radial
cells and 12 multipoles, moving the matching rectangle from 2 to 4 m changes
source work by **0.00206%**. The corresponding fixed-wall change is **23.646%**,
using the 4 m result as denominator for both percentages. The prior roughly 19%
figure used the 2 m fixed-wall result as denominator. Open pulse energy residuals
are at most 2.21e-14 relative to initial response energy; mass change is at most
2.83e-16 and modal stationarity residual at most 6.32e-13. Beyond-grid magnetic
energy is positive and included. No exterior fluid or state repair is introduced.

The pulse still is not spatially converged: open-boundary final peak speeds are
0.05919, 0.13227, 0.16262 and 0.18119 m/s at 8/16/24/32 cells. At 16 cells,
halving CFL changes source work by about 0.046%, larger than the measured boundary
surface effect. More modes alone cannot fix fluid or time-discretization error.
The current single-center expansion requires a matching rectangle enclosing the
entire chamber's source sphere; efficiency for long chambers remains to evaluate.

Derivation, measured results, inspected figure, complete inputs/histories and
source fingerprint are in `docs/OPEN_BOUNDARY.md`,
`docs/OPEN_BOUNDARY_VERIFICATION.md`, and
`outputs/open_boundary_verification/metrics.json`. This verifies the open
**magnetic** exterior, not a moving free plasma surface, gas outflow, transport,
magnetic shocks or FRC formation. `frc-vacuum` and the full objective stay active.

## Magnetic shock and implicit solver checkpoint

That checkpoint had **99 passing tests** (25.17 s). That continuation made
verified progress; all measurement/report jobs completed, with none left live.

An independent strong perpendicular shock now has exact Rankine–Hugoniot jumps
and steady cylindrical branches. Finite-difference momentum and invariant checks
validate the reference independently of the numerical kernels. Both actual 2-D
solvers are tested in an interior annulus before remote boundaries dominate.
The centered coupled method creates a **false Bz reversal** at nr=48, reaching
Bz/B1=-0.027288 even with a relative response-energy residual around 1e-14.
Such magnetic structure must not be interpreted as FRC topology.

The explicit `shock_diffusion` option adds a nonnegative compressive diffusivity
that scales as h² in smooth flow and h at shocks. It is off by default and
separately records numerical magnetic heat, without changing physical Ohmic
heat or putting fluid/resistivity in vacuum. No cell floors or energy repair are
introduced. With C=1, magnetic-profile L1/B1 errors at nr=24/48/96/192 are
0.330402, 0.218542, 0.106434, 0.0581267; density errors are 0.299483, 0.187281,
0.0960794, 0.046949 kg/m³; pressure errors are 1.93626, 1.06167, 0.574753,
0.190947 Pa. No false sign reversal occurs in those C=1 windows. However,
frozen B/rho error is **not monotonic**, and nr=192 still overshoots Bz/B1 to
3.6293 versus the physical shock value 3.53234. The conservative bounded
reference preserves that special B/rho invariant to roundoff. More diffusion
is not uniformly more accurate; the coefficient study retains that tradeoff.

At nr=48, moving the inner nonanalytic extension changes pressure/B profiles by
0.3302%/0.08888% in relative L2. Doubling chamber length changes profiles below
7.1e-9; a shorter 0.01 s wall-taper control gives density change 4.86e-5.
CFL halving changes the B profile by 1.43e-5, then 3.45e-6. These are measured
controls at that grid/window/time, not general free-boundary verification.
Maximum coupled response-energy residual is 1.94e-14, mass change 2.19e-16,
and face divergence 1.25e-16 T/m. Numerical magnetic heat remains explicitly
reported (0.433 J at nr=192 for the whole chamber); physical eta and heat are
zero. Entropy-producing shock dissipation need not vanish with grid spacing.

The midpoint solve now defaults to a cached-vacuum Picard iteration at tolerance
1e-13. It solves the same implicit equations and retains atomic failure. The
original direct solve remains selectable. Direct/Picard tests cover both
magnetic boundaries and independent physical/numerical diffusion choices. Their
nr=24 shock profiles agree within 9.9e-13 relative L2; measured advance time is
6.96 times faster for Picard in that fixture. This is a CPU wall-time measurement,
not a GPU or general real-time claim.

The inspected six-panel figure, 20 completed cases, complete configuration,
initial/final diagnostics, profiles, environment and source fingerprints are in
`docs/SHOCK_VERIFICATION.md`, `docs/SHOCK_NUMERICS.md`, and
`outputs/shock_verification/metrics.json`. Old reports retain historical hashes.
Next, investigate the split method's inconsistent magnetic/fluid numerical flux
and local momentum/weak-shock accuracy, while preserving the verified vacuum
and work identities. Continue the unresolved cold rarefaction, equilibrium,
skin-layer, free-interface, transport and formation requirements below.

## Conservative fluid/CT coupling checkpoint

The current full suite has **109 passing tests** (32.92 s). This continuation is
verified progress: all 22 measurement jobs completed, including 20 successful
runs and two deliberately disabled-method failure controls. No job remains live.

`ConservativeVacuumMHD` advances mass, momentum, total fluid-cell MHD energy and
upwind constrained induction together, with shared Rusanov/UCT diffusion. It
retains the exact discrete induced vacuum and infinite magnetic exterior, with
no exterior fluid or conductivity. Vacuum-side magnetic traction and Poynting
flux include interface-current work. Positive numerical diffusivity D=a h/4 is
assigned only to plasma boundary nodes and uses the symmetric stiffness current;
its diagnosed heat is separate from physical eta J². This is declared numerical
interface stabilization, not transport. Physical and numerical electric fields
enter all adjacent total-energy fluxes, including transport along the boundary
layer, with no double-counted Joule source. An independent manufactured Poynting
identity checks that a prescribed electric curl in initially uniform B changes
magnetic energy without a spurious instantaneous thermal-energy derivative.

The original strong shock fixture reaches 0.02 s at nr=24/48/96/192. Interior
frozen B/rho defects are 2.85e-6, 9.22e-8, 4.45e-13 and 9.27e-14, with no false
field reversal. B-profile L1/B1 errors are 0.309278, 0.203389, 0.104325 and
0.0746278. Cold boundary cells require blending complete conservative/CT fluxes
toward an admissible first-order Euler proposal: 0/88, 0/176, 6/372 and 154/810
Euler stages are limited. Minimum high-order weight is 0.000332 at nr=192.
No cell pressure/density floor or energy repair is used. A failed low-order
proposal still raises atomically; no unrestricted positivity theorem is claimed.
The disabled limiter fails at 0.019111 s on nr=96; disabling interface diffusion
and limiting fails at 0.0100 s on nr=24. These failures remain in the data.

The invariant improvement is not universal profile improvement: at nr=192 the
new B error 0.07463 exceeds the historical split C=1 error 0.05813, and pressure
error is 0.36880 Pa versus 0.19095 Pa. Global limiting caused by boundary cells
reduces accuracy in the whole domain, so localization is the immediate numerical
task. Shock time halving at nr=96 changes B/density profiles by 0.02976% and
pressure by 0.07051%; doubling chamber length at nr=48 changes profiles below
2.70e-8. The original initial data and independent reference were retained.

Fluid-cell energy flux accounting is at roundoff, but independently calculated
gas-plus-full-space-induced-field response energy is not forced exact. Relative
shock response errors decrease 1.996e-4, 8.248e-5, 4.050e-5, 2.873e-5. This
explicit finite-grid error is distinct from the historical split work identity.
Mass changes are at most 2.19e-16 and compatible divergence 2.09e-16 T/m across
completed cases. All quantities remain independently diagnosed.

The same prescribed-coil pulse runs at nr=8/16/32/64. Final response error
normalized to **actual source work** decreases in magnitude 2.722%, 1.108%,
0.3406%, 0.1021%. At nr=64, physical Ohmic heat is 0.082397 J and numerical
interface heat is 0.0014725 J (about 1.8% of physical heat). Source work still
changes 0.125487→0.121817 J from nr=32→64, and peak speed increases
0.06981→0.12737→0.17880→0.20744 m/s over the four grids. The driven flow remains
unconverged. No pulse stages require flux limiting. CFL halving at nr=32 changes
work by 0.001741%, and magnetic matching extent 2→4 at nr=16 changes work by
0.002050%. Removing numerical interface diffusion at nr=32 changes work by
0.6537% and peak speed by 1.4146%, so its effects must not be ignored.

Independent smooth advection errors at nr=12/24/48 are 0.05585, 0.009702 and
0.001778; Gaussian diffusion errors are 0.007897, 0.002115 and 0.0004833, with
no flux limiting. Harmonic-polynomial tests independently verify vacuum-side
interface traces and their radial convergence. Tests cover no vacuum numerical
resistance, rest/bias preservation, cold-boundary atomic failure, local energy
conservation and independent source-budget refinement. The inspected figure,
full inputs/histories/profiles/failure records and fingerprints are in
`docs/CONSERVATIVE_COUPLING.md`, `docs/CONSERVATIVE_VERIFICATION.md` and
`outputs/conservative_verification/metrics.json`.

This remains a fixed gas-filled chamber reference. Free plasma motion and
actual end conditions, transport/ionization/Hall validity, converged formation
and a dynamic experiment UI remain unfinished. Original Stage 2 stays gated.

## Local complete-update limiter checkpoint

This continuation is verified numerical progress, not completion of the FRC
objective. The suite has **117 passing tests** (46.65 s). The default conservative
solver now chooses first-order shared face/corner fluxes only around inadmissible
Euler cells and rechecks newly troubled neighbors. Metric sources and actual
boundary energy fluxes follow that same local choice. Every flagged cell has
the independently checked low-order update, giving conditional finite
termination without a hidden global fallback. No floors or thermal repairs
were introduced; failed baseline proposals still raise atomically. The previous
global convex limiter remains an explicit, freshly reproduced control.

On the unchanged original shock, nr=192 B-profile L1 error decreases from
0.0746278 globally to 0.0593295 locally, and pressure L1 from 0.368802 to
0.334502 Pa. Frozen B/rho defect remains 8.06e-14. Only the two outermost radial
cell layers are directly affected; maximum stage support is 0.8545% of cells,
with zero direct replacements in the central independent-reference window.
Local limiting occurs more often (320 stages) but leaves distant fluxes intact.
The manufactured neighbor-chain test requires four successive detection passes;
the original nr=192 case needs at most two. Mass, fluid energy-flux accounting,
compatible divergence and vacuum stationarity remain near roundoff.

This improvement exposes a separate accuracy limitation. The nr=192 local
response-energy residual is +0.0140806 J (1.40031e-4 of initial response energy)
versus -0.00288898 J globally. The local residual is also **5.22% of the net
outward fluid energy**, a more demanding transfer scale, though that boundary
transfer is not external source work. Halving CFL raises the local residual to
0.0237089 J while central B/density profiles change only 0.001610% and pressure
0.002069%. Thus the interface evolution/energy is not time-converged simply
because the central shock is. Doubling chamber half-length at nr=96 changes
central profiles below 9.01e-9; the independent response budget remains separate.

The additional nr=384 run completes to the original 0.02 s end time: central
B/density L1 decreases to 0.0306517, pressure L1 to 0.175216 Pa, frozen-ratio
defect is 9.72e-13, and no reference-window cell is directly limited. Maximum
affected support is 0.4286% of cells, with at most two replacement passes.
However, response-energy error grows to **0.0575258 J, 30.20% of net outward
fluid energy**, and numerical interface heat rises from 0.06331 to 0.09848 J
between nr=192 and 384. The interface energy budget is therefore not converged.
Minimum density reaches 0.000388047 kg/m³ without a floor; 1575 accepted steps
take 2064 s after 101 s setup on the recorded CPU environment. The large-grid
run finished; no numerical process remains pending at this checkpoint.
All 21 case records (20 completed, one deliberately disabled-limiter failure)
and both energy audits have been source-fingerprint checked, and the final
scientific figure has been visually inspected.

The same coil pulses, time/matching controls, and smooth advection/diffusion
cases reproduce historical measured source work, response energy and analytic
errors exactly, with no local limiting. The earlier 0.1021% pulse error relative
to source work and unresolved peak-flow/skin-layer limitations remain.

An instantaneous response-rate identity now separates outward fluid power,
cell magnetic-energy rate and full-space nodal magnetic work. An independent
symmetric perturbation completes the magnetic vacuum and differences response
energy. It checks this diagnostic without changing any evolved state or budget.
The best symmetric-difference errors across the five snapshots stay below
2.62e-8 W at nr=96 and 3.65e-7 W at nr=192; probing leaves final diagnostics
exactly unchanged. At the final nr=192 state, the high-order defect rate is
-204.075 W, the first-order rate -2.83969 W, and the selected local rate
-7.23180 W. Thus this wall-state discrepancy is present in the spatial update,
well above diagnostic roundoff; its cumulative and timestep-sensitive behavior
still needs a consistent interface treatment.
The derivation is in `docs/INTERFACE_ENERGY_AUDIT.md`; measured rates, perturbation
errors, full controls, limiter maps and source fingerprints are retained in
`docs/LOCAL_LIMITER_VERIFICATION.md` and `outputs/local_limiter_verification/`.
The next interface task is to use that decomposition to improve consistent
vacuum/CT energy exchange without sacrificing local conservation or adding heat
repairs. Stronger multidimensional shocks, cold rarefaction, free surfaces,
actual ends, transport/ionization/Hall validity, formation and the dynamic
experiment UI all remain in scope.

## Compatible-energy and free-boundary-reference checkpoint

Positive staggered face-square magnetic energy and a natural variational
vacuum trace are now explicit solver options, with optional SSP RK3. Pressure,
initial energy, complete-update admissibility and diagnostics use the same
selected quadrature. Boundary work uses actual selected CT corner fields;
there are no thermal repairs or exterior fluid variables. Historical defaults
remain selectable and unchanged. Independent polynomial/force/work tests cover
unequal grid spacings, nonautonomous temporal order, third-stage atomic rejection
and mixed-boundary fallback. The pulse16 and advection12 legacy checks reproduce
historical source work, energy and analytic errors exactly.

The original shock completes nr=24/48/96/192/384. At nr=192 the response residual
falls from historical 0.0140806 J to 0.000596124 J; at nr=384 it falls from
0.0575258 J to **0.00173278 J**, or 0.728851% of net outward fluid transfer.
The fine B L1/B1 remains 0.0306581 and pressure L1 0.175312 Pa; frozen-ratio
defect is 7.07e-13. No central-reference-window cell is directly limited.
The finest run finishes 1408 steps with minimum density 0.000916625 kg/m³,
3423/4224 limited Euler stages, at most 0.48489% affected cells per stage and
three replacement passes. It takes 98.9 s setup and 2792.9 s advancement.
Its energy error still grows from nr=192, so full boundary convergence is
not established by this substantial improvement.

A per-step quadratic RK audit independently predicts the entire nr=192
response error to roundoff, with maximum stage spatial defect 2.20e-11 W and
maximum per-step identity error 1.13e-12 J. The constant-bias semidiscrete work
mismatch is resolved for this method; its finite-step error is not repaired.
At nr=24, CFL halving reduces RK3 energy error by factors 8.88 and 7.85.
The fine boundary is not time-converged: half-CFL gives 0.00100662 J with RK3
and 0.00218710 J with RK2. The full first-order spatial baseline instead gives
3.97637e-7 and 4.79246e-8 J, but degrades B L1 to 0.15656.

An explicitly experimental fixed first-order boundary strip retains the
second-order interior. Recorded coarse/fine controls test that hypothesis;
its forced order reduction is counted separately from any adaptive fallback.
It is not a production default. Source-driven skin layers, boundary accuracy
and finer shock meshes must be checked before adopting it.
The nr=64 unmodified compatible pulse residual is **0.00315436%** of applied
source work, versus historical 0.102119%. Peak flow remains mesh-sensitive.
All original pulse/smooth cases complete with inactive adaptive limiting.

The spherical FRC reference independently verifies true vacuum outside,
continuous magnetic traction with no sheet current, exact force balance and
full-space energy integrals. The existing topology detector finds one O-ring
and two axis nulls; separatrix-volume error falls from 0.219889% to 0.00374111%
over 24→192 intervals. This is initialized ideal equilibrium, not formation,
a density/transport closure, stability or confinement evidence. Its figure
has been visually inspected. See FREE_BOUNDARY_REFERENCE.md.

The full suite now has **135 passing tests (49.04 s)**. All 25 primary cases
(24 complete, one retained disabled-limiter failure), two first-order controls,
six boundary-strip probes and four RK audits are finished. Their source hashes
match the frozen engine/measurement fingerprint
`5247cc3d60224324a8f1f731b164f9fd41a16a8ea70a140af70ae1a302c022c2`.
Both scientific figures were visually inspected. No numerical jobs remain.
The nr=192 strip probe preserves B L1 0.059316 and gives response residual
3.97369e-7 J, decreasing to 4.78879e-8 J at half CFL; its 766 forced cells
and 1528 affected cells are explicitly reported separately from adaptive
fallback. Source-driven and finer-boundary checks remain open for this probe.
See COMPATIBLE_ENERGY.md and COMPATIBLE_ENERGY_VERIFICATION.md; full
source-matched records are in outputs/compatible_energy_verification.
The next numerical step is to verify the fixed boundary strip with driven
pulses and finer meshes, or improve the higher-order boundary closure without
losing conservation, CT and admissibility. Free-boundary evolution, actual
ends, transport/ionization/Hall validity, converged formation and the dynamic
experiment UI remain required. The Stage 2 formation gate stays closed.

## Boundary source, skin and end validation checkpoint

The current suite has **142 passing tests (57.97 s)**. A new independent
infinite-cylinder Bessel solution checks sinusoidal skin diffusion, regular-axis
and vacuum matching, Maxwell/Ohm identities and time-averaged surface/volume
Ohmic power. Fifteen recorded grid/time cases compare physical-only diffusion,
the production induction stencil and the fixed first-order boundary strip.
Both production choices converge in magnetic field; the surface-current
stencil has separate, slower pointwise accuracy. At nr=128 the maximum B L2
errors normalized by applied amplitude are 0.000934414 (standard) and
0.00150768 (strip), with about 7.22 cells per analytic skin depth.
This is stationary-conductor induction verification, not evolved gas dynamics.

All twelve full coil-pulse runs and four gas-end controls are complete, with
hashed lossless conserved/primitive/CT/current snapshots at 0.04, 0.06 and
0.10 s. The standard nr=16/32/64 budgets and peak speeds reproduce the previous
compatible study exactly. At nr=128, standard/strip final peak speeds are
0.221845/0.217821 m/s, but their maximum same-grid velocity L2 difference is
**9.76713%**. The 64→128 velocity changes are **7.05911% / 10.4227%**.
Thus the strip's favorable shock energy result has a material driven-flow
accuracy cost and it remains experimental; no production default is changed.
Density/pressure perturbations and surface currents also remain resolution
sensitive. A small energy residual is not an accurate formation result.

Half timestep changes pulse velocity below 0.0000710%; moving the magnetic
matching surface changes it below 0.000333%. Moving actual reflecting gas ends
from ±1 to ±2 m at a common matching extent of 3 m changes central |z|<0.5 m
velocity by 0.8174% / 0.8372%, and whole-domain source work by about 2.1%.
Spatial and material-boundary accuracy dominate those time/matching errors.
The change of gas ends is a physical-domain sensitivity, not an open-end model.

The nr=384 original shock with the strip completes to 0.02 s: B L1/B1 is
0.03065827, pressure L1 0.1753145 Pa, frozen-ratio defect 5.11e-13 and response
residual **1.01130e-7 J**, versus 0.00173278 J for the standard compatible case.
There are 985 steps and no additional adaptive limiting, but 1534 forced and
3064 neighboring affected cells. Minimum density is 0.0496321 kg/m³ without
a floor; fluid energy-flux residual is 3.28e-15 J. This preserves interior shock
accuracy and removes the observed growth of the strip's energy error through
nr=384. A separate nr=384 time study is not claimed.

All **32 new recorded cases** (12 pulses, 15 skin cases, four end controls and
the fine shock) are complete. No numerical jobs remain. Full fields, inputs,
source/environment fingerprints and the report are in
`outputs/boundary_verification` and `docs/BOUNDARY_VERIFICATION.md`.
An independent annular polynomial checks field restriction to 4.44e-16;
all saved-field hashes and the common frozen engine fingerprint match.
The final six-panel scientific figure was visually inspected after fixing
overlapping resolution-axis labels.

The next numerical implementation must address boundary/flow accuracy and
free-interface evolution using these source/skin/current checks. Keep the strip
experimental until its driven-field accuracy is adequate. Real gas ends,
transport/ionization/Hall validity, stronger shock/cold-flow and curved-equilibrium
accuracy, converged formation and a dynamic experimental UI remain in scope.
The Stage 2 formation gate and full objective remain unchanged.

## Moving material / infinite vacuum checkpoint

That checkpoint had **155 passing tests (61.39 s)**. That continuation was
verified progress toward the active objective, not completion. All 26 evolution
cases and 27 static evaluations completed. No numerical jobs remain live.
Source/configuration/environment and all lossless field hashes are checked;
both scientific figures were visually inspected, with logarithmic tick-label
collisions corrected. The full report is MOVING_PLASMA_VERIFICATION.md and the
derivation is MOVING_PLASMA.md, with records in outputs/moving_verification.

The new `MovingPlasma` reference has actual two-dimensional material motion,
adiabatic gas only inside its evolving surface, regular a=psi/r² induction
labels, analytic magnetic/pressure virtual work, an exactly curved matching
sphere and decaying magnetic modes extending to infinity. No exterior fluid,
pressure/density floors, numerical dissipation or energy repair are introduced.
This prototype is ideal-only with prescribed uniform bias; fixed-chamber
resistive references and production defaults remain unchanged.

Independent uniform-sphere/dipole energy and homologous force errors converge
approximately at second order. At nr=32 and scale=0.9, relative energy/force
errors are 0.000724094 / 0.002616710. Translating the sphere by 0.3 m exposes
spurious net axial force -0.00796587→-0.000122393 N at nr=4→32; it refines
but precludes a claim of exact magnetic axial-momentum conservation.

For exact free gas expansion at 0.1 s, velocity relative mass L2 errors are
0.0726360, 0.0268221, 0.00924936 and 0.00302204 at nr=4/8/16/32. Finest
surface-position error is 0.000138224 m and pressure-cell-average relative L2
error 0.000101608. Actual material mass is unchanged. The exact sphere radius
is 1.009950494 m; polygonal geometry and cell entropy projection are not repaired.

For the initialized spherical FRC, initial acceleration mass RMS decreases
0.713087→0.243325→0.0828970→0.0285246 m/s². Maximum artificial speed after
0.04 s decreases 0.0764241→0.0327842→0.0120209→0.00635480 m/s; finest surface
displacement is 0.000127825 m. This is measured equilibrium drift, not formation
or stability. Moving the matching sphere 2→3 m at nr=16 changes the artificial
velocity by at most 0.000331791 m/s RMS (12.2083% relative to the small reference
drift). Quadrature 5→7 changes it below 1.65e-12 m/s RMS; modes 4→8 below
2.15e-6 m/s. A tiny vacuum residual does not erase matching/mesh error.

Velocity Verlet and source-work integration show second-order time refinement
when the explicit cap controls actual steps. In the Gaussian uniform-bias ramp,
nr=8 response errors at dt=0.001/0.0005/0.00025 s are 0.000448658,
0.000112158 and 0.0000280390 J. However driven velocity changes **47.1666%**
between nr=8 and nr=16, and peak speed rises 0.166450→0.290958→0.524906 m/s
at nr=4/8/16. The driven ideal current-layer response is explicitly unconverged.
Source-work accuracy is not a reason to interpret it as physical formation.

Vertex psi is preserved near roundoff, while interpolated interior material
flux is not exact. A prescribed nonuniform deformation has relative flux errors
0.00386486→0.000929701→0.000229928→0.0000573222 at nr=4/8/16/32. Initialized
FRC interior flux drift at nr=32 is 1.85284e-6; source nr=16 is 0.000631995.
These are numerical interpolation effects, not physical resistivity or
reconnection. The regular-axis stretch approximation and low-energy/distorted
mesh modes remain explicit. No shock capture or remapping is provided yet.

The next resistive checkpoint below implements the scalar diffusion and heat
extension against the independent spherical-decay eigenmode. This ideal report
retains its historical source, including ideal source-flow limitations.

## Resistive moving material / vacuum checkpoint

The full suite has **170 passing tests (62.11 s)**. All 13 fixed-geometry decay
cases, independent stiff-mode controls and 20 actual moving cases completed.
Full configurations, histories, snapshots, source hashes and inspected figures
are in `outputs/resistive_moving_verification`; the method and report are
`docs/RESISTIVE_MOVING_NUMERICS.md` and
`docs/RESISTIVE_MOVING_VERIFICATION.md`. Frozen production/measurement SHA256 is
`452e51b4405b8a44fc6191f6716fa4e19e5798ba367bb84a0953ce3c6703503b`.
Report, analysis and plotting sources have separate recorded hashes.

`MaterialResistivity` uses a consistent plasma-only weighted mass and the same
full-space magnetic stiffness to solve scalar physical induction and recover
Jtheta. Implicit midpoint gives positive per-cell eta J² heat and a separate
source-work/energy check. `ResistiveMovingPlasma` atomically composes two
resistive half-steps with ideal material motion, updating both flux labels and
cell entropy. Fixed mass is retained and no exterior fluid, artificial
resistivity, pressure floor or energy repair is introduced. Eta=0 follows the
unchanged ideal evolution. Failure of the second resistive stage after trial
heating/motion leaves all accepted arrays, time and accounting unchanged.

An independent continuum conducting-sphere eigenmode matches a vacuum dipole
with continuous field at the surface. Maxwell/Faraday/Ohm derivatives and
all-space quadrature independently verify its current, energy and heat. At
t=0.1 s and nr=4/8/16/32/64, magnetic-field relative volume L2 errors decrease
0.0926223→0.00523351, current errors 0.0521209→0.000206685, and local-heat
errors 0.124056→0.000458581. Q1-gradient field error is about first order;
flux, current and cell heat are about second order. Axis-field and one-sided
surface-current errors are reported separately. A second radial eigenmode,
time, tolerance, matching-radius, mode-count and quadrature controls are complete.
The geometry is held fixed for this diffusion fixture; it is not force balanced.

A raw midpoint step of an independently diagonalized decaying spatial mode
produces a false amplitude -1/3 instead of exp(-4), despite a 2.22e-16 J energy
residual. Adaptive step doubling controls field, current and local heat, accepts
the unmodified two-half solution, and retains positive heat. Tolerances 1e-3
through 1e-6 reduce absolute modal-amplitude error 3.30e-4→4.13e-6 and preserve
the tested positive decay. This is local error control, not a global error,
positivity or arbitrarily small-field topology guarantee.

The actual moving cases start from the analytic FRC and then undergo resistive
relaxation, with or without a Gaussian uniform-bias ramp. Reversed topology is
initialized, not formed. With D=eta/mu0=0.1 m²/s, the driven nr=16→32 comparison
differs by 63.821% in mass-weighted velocity, 29.357% in magnetic change,
3.67623% in total current and 5.4924% in local heat density. Density change
differs by 54.0889%. These material-coordinate differences use coincident
material sampling and actual volume restriction, independently checked with
polynomial and geometric integration tests. They are not absolute solution
errors. Higher physical D=1 still gives 53.915% velocity difference; changing D
is a physical change, not numerical stabilization. The large flow error remains.

At nr=8, halving outer time caps gives approximately fourfold reduction in
velocity differences and response residuals. Driven energy residuals at
dt=0.001/0.0005/0.00025 s are 4.39132e-4/1.09776e-4/2.74435e-5 J; the last
velocity comparison differs by 0.0149697%. At nr=16, moving the magnetic
matching radius 2→3 m changes driven velocity by 0.539718%; mode/quadrature
changes are smaller. Tightening diffusion tolerance 1e-4→1e-5 leaves identical
accepted stages in these D=0.1 cases and therefore identical results. D=1 cases
do reject/subdivide diffusion trials. Mass is unchanged; maximum sampled
resistive-stage energy residual is 1.679e-14 J and vacuum relative residual
4.451e-13. These are distinct from full moving-state accuracy.

Next inspect spatial modes and surface thermal/force accuracy and add circular
coils at moving coordinates with the geometric external-source derivative.
An independent exact dipole/loop force fixture is specified in MOVING_COIL_PLAN.md.
Shock capture/remapping, material distortion, physical transport/ionization/Hall
regime, actual ends, converged formation and the experimental UI remain required.
All overarching issues and the Stage 2 entry remain open.

## Surface-resolution and thermal-compatibility checkpoint

This continuation implements optional `radial_grading` in `MovingDomain` and
records actual angular resolution in the resistive case contract. The default
uniform implementation is preserved: all archived arrays in the new source
nr=8/16/32 baselines are bitwise identical to the preceding checkpoint. The
frozen production/shared-measurement fingerprint is
`2749bf18bff1bdab93b1c338b7944ec7a6773d23c39076700f83a4faec9b90c8`.
All 21 dynamic and 13 independent cases completed, including the supplemental
finest quarter-time control (798 accepted motion steps). Analysis/render and
supplemental reproduction helpers are separately hashed. Full configurations,
lossless snapshots, histories, comparisons and figures are retained in
`outputs/surface_verification`; see SURFACE_NUMERICS and SURFACE_VERIFICATION.

Nine new independent tests bring the suite to **179 passing tests (80.30 s)**.
Revolved-polygon moments check graded geometry; exact free-gas expansion and
conducting-sphere decay check the actual kernels. Current/local heat remain
approximately second order under refinement. Grade-2 nr=32 current and local
heat relative errors are 0.000829269 and 0.00182061; Q1-gradient B error is
0.0105070. Grading trades interior resolution for surface resolution and does
not automatically improve the smooth free-gas reference.

Material-field sampling now inverts the actual initial sector geometry instead
of comparing equal normalized coordinates on different polygons. Exact convex
cell intersections resolve discontinuous thermal fields before positive Gaussian
integration. Independent affine field and geometric/indicator integrals verify
this diagnostic. Across all 18 comparisons, increasing integration order 4→6
changes reported relative norms by at most 3.38e-9 of the norm. This removes
the observed diagnostic-quadrature uncertainty without smoothing evolved fields.
The velocity norm integrates actual Q1 fields with initial physical mass;
thermal fields are compared pointwise rather than conservatively averaged onto
coarse cells. These definitions differ from the preceding report and are explicit.

The new velocity differences are 34.1273% for uniform nr=16/32 and 55.2917%
for grade-2 nr=16/32. Fine grade-2 angular doubling changes velocity 15.8362%.
Successive fine time controls change it 5.39874%, then 2.75800%; the ratio near
1.96 and nonmonotonic current differences do not establish asymptotic order.
Boundary-speed differences remain 0.04457785 then 0.01810780 m/s. The finest
quarter-time run has peak speed 0.5668651 m/s, minimum density 0.0555440 kg/m³,
physical heat 0.2738125 J and response residual 1.97602e-6 J. No floor, energy
repair or change of physical D is used. Spatial and time uncertainty remain.

Saved-state force decomposition exposes cancellation between large pressure
and magnetic surface accelerations. Recovered weak-current J cross B is a
separate discrete diagnostic, not an exact continuum reference. Mass remains
unchanged; maximum sampled resistive-stage energy and vacuum residuals are
1.386e-13 J and 7.566e-13. These do not establish local boundary accuracy.

The initial thermal-compatibility audit independently recovers the boundary
current from one-sided Maxwell derivatives with relative error below 4.7e-11.
The pressure equation gives positive equatorial initial p-dot of 3.75 Pa/s for
the unforced FRC and 1.379547904 Pa/s for the source case; a persistent smooth
p=0 trace under continuous-B/no-extra-stress assumptions would require zero.
The exact conditions and primary sources are in SURFACE_COMPATIBILITY. Keep
these initialized cases as stress tests, investigate a justified edge closure
or verified nonsmooth expansion, and do not substitute arbitrary eta tapers,
removed edge heat or easier fixtures for physical formation validation.

Both final scientific figures were visually inspected; colliding logarithmic
minor x-axis labels were corrected. The strict report serialization was fixed
to write diagnostic arrays as numeric lists, then the complete assembly was
rerun using all 18 verified comparison caches. All jobs are terminal. The goal
and relevant beads issues remain open with physical surface closure and the
independent spatial-coil force reference as concrete next work.

## Electron-current inertia and conditional transport checkpoint

This continuation adds `MaterialCurrent` with physical me and frozen prescribed
ne/eta profiles, counting current kinetic energy alongside magnetic energy and
positive per-cell collisional heating. Its sparse mixed Ampere/Ohm solve keeps
vacuum elliptic with no exterior carriers. Increment formulation preserves
stationary current; diagonal algebraic scaling improves axis-conditioned
roundoff without changing physics. A final exact-ideal-limit regression also
preserves tiny material flux under a much larger source-field change, with its
Ampere residual explicitly checked. Independent massless-limit and variable-
coefficient physical-volume energy checks pass. This is a separate linear
component, not a hidden change to the moving plasma or a complete edge model.

Thirteen new tests bring the full suite to **192 passed (79.10 s)**. Primary NRL
collision coefficients in SI agree with independent published cgs numbers.
Coulomb logs are never floored; invalid classical closure is explicitly
unavailable. Reference Stage 0/1 metadata now carries conditional thermal,
Hall/electron-inertia/Debye and drift diagnostics, retaining the original
prescribed conductivity and distinguishing inferred Drude rates. No absence of
a warning establishes validity; neutrals/ionization and tensors remain open.

All 35 current cases completed, including the added thin-layer nr=64 control.
Uniform inertial sphere decay, collisionless modified-Bessel screening and
an exact variable-density/resistivity profile independently test continuum
fields/current/heat. For de/R=.2, nr=8/16/32/64 current errors are
0.0249752/0.00768733/0.00247993/0.000828866; finest B and local-heat errors are
0.00523256 and 0.00104859. Finite inertia retains the weak derivative error in
the interpolated initial field, yielding roughly order 1.5 current convergence;
surface grading alone worsens this in the recorded control. This is explicit.

Strong screening de/R=.05 exposes false reversal on nr=8: minimum Bz/bias
-0.00330320 and 40.2668% of sampled material volume negative, while the analytic
field is positive everywhere inside. The amplitude is small relative to bias
but would invalidate a sign-only topology claim. Finer recorded samples have no
negative Bz; nr=64 current error is 0.00480722 and magnetic-field error 0.0730778.
No clipping, carrier floor, anomalous resistance or heat cutoff removes the
coarse artifact. Sample positivity is not a between-sample topology guarantee.

An independently diagonalized fixed spatial mode gives temporal halving ratios
4.002221 and 4.000555. Matching-radius, angular, quadrature, mode-count and
grading controls are retained. True-sphere energy quadrature orders 48/96 agree
within 9.172e-14 relative; maximum whole-run energy residual is 8.414e-15 of
the reported continuum energy scale and maximum equation residual 2.341e-12.
These do not erase current/field error or validate nonlinear formation.

The conditional scan records 43 assumed states and 41 prescribed-density-edge
samples. The original Ar+ case gives di/R=4.09721, de/R=0.0151831 and electron
mean-free-path/R=0.110877. Conditional perpendicular momentum and Z=1 parallel
resistivities are 8.57541e-5 and 4.37521e-5 ohm m; their limits are distinct and
neither is automatically installed. Fixed charge/temperature/ionization are
inputs, not equilibrium predictions. The prescribed vanishing-density edge
and its fixed current are a diagnostic illustration, not a thermal solution.

All source/configuration/lossless archive identities and companion documents
are checked. Frozen component measurement SHA256 is
`c66fb837ef6af44768c87d3a88e624f9ec75e318fd24a4cfd62a82752725908f`;
report, plotting, time audit and transport scan are separately hashed. Both
final figures are visually inspected. Full records are in
`outputs/current_verification`, equations in CURRENT_INERTIA and
TRANSPORT_REGIME, and results in CURRENT_VERIFICATION. All jobs are terminal.
EXTENDED_PLASMA_PLAN defines the next Hall/toroidal wave reference and the
remaining current-energy, thermal and edge integration requirements. The full
goal remains active; a component milestone is not FRC formation validation.

## Linear Hall / toroidal reference checkpoint

The full suite now has **208 passing tests (78.87 s)**, including 16 independent
Hall/finite-mass tests. `hall_parameters`, `hall_modes` and `hall_cylinder` add a
bounded uniform-density incompressible reference with regular poloidal/toroidal
fields and velocities. Exact center-of-mass/reduced-mass coefficients follow
parent species momentum; drift energy, bulk/magnetic exchange and positive heat
are paired. There is no exterior gas and no claim that its restricted periodic
reference boundary is a real plasma/vacuum or perfect-conductor boundary.

All 49 cylinder cases and 16 parallel modes (three time steps each) completed.
Actual cylinder field errors decrease from approximately 5.2% at nr=8 to 0.626%
at nr=64. Both helicities and prescribed scalar resistive damping are measured;
physical finite mass recovers ion/electron cyclotron asymptotes in the parallel
limit. A poloidal perturbation generates toroidal field and flow through Hall;
the MHD control does not. Physical volume integrals independently verify all
three energies, weak current and local midpoint heat; angular momentum is
conserved in this reference. The original argon target has not been replaced
by the separately labeled proton verification fixtures.

Parallel time-halving ratios are 3.998876–3.999721. Initial interpolation in the
slow cylindrical branch excites small fast secondary modes, giving nonmonotonic
apparent projection damping at coarse time steps. An independent full-spatial
helicity exponential and additional 2048/4096/8192-step controls resolve this:
relative time errors decrease 1.58885e-7 → 3.97212e-8 → 9.93030e-9. Remaining
spatial/projection error is explicitly separate. A one-step collisionless mode
still has 21.46% phase error despite roundoff energy conservation; no energy-only
accuracy claim is made.

Maximum whole-run energy-plus-heat residual is 6.63520e-15 relative, equation
backward error 2.38842e-17 and independent physical/modal energy disagreement
5.12417e-14 of initial energy. Cylinder source fingerprint is
`868ece2de1730a51d4aaf363652e3a19139c8289e046eda3188b57a0a45b0e3f`.
Full source/config/archive/history checks pass; records and inspected scientific
figures are in `outputs/hall_verification`. Equations and findings are in
LINEAR_HALL_NUMERICS and LINEAR_HALL_VERIFICATION. All jobs are terminal.

Next implement TOROIDAL_INTEGRATION_PLAN's conservative angular momentum and
toroidal induction in the bounded compressible core, verify geometric rotation/
pinch forces and nonlinear transport, then integrate Hall work and finite-mass
thermal/current transport from the parent species equations. The current
linear reference must not be inserted as an incomplete nonlinear closure.
Heated edge physics, spatial coils, actual vacuum/ends, shocks/remapping,
kinetic ordering, converged formation and experimental UI remain required.
The full goal remains active, and Stage 2 remains guarded.

## Nonlinear compressible toroidal checkpoint

The full suite passes **224 tests (85.76 s)**. Six new production modules add a
separate compressible state/solver with toroidal field, rotation, conservative
axial angular momentum and full magnetic/centrifugal stresses. Exact annular
r² moments enter pressure recovery and rotation/field energy, including the
axis. Both first/second-order zero-toroidal RHS recover the existing poloidal
solver, whose API and historical results remain unchanged. Failed final stages
leave accepted state, time, heat and limiter accounting untouched.

Scalar toroidal resistivity uses a positive compatible tridiagonal current mass
and matched Faraday/Poynting/heat work. It exactly recovers the independently
known 8Db diffusion of g=a+b r² at the axis and wall. Axis Ez is solved rather
than incorrectly set to zero; the conducting outer wall still has Ez=0. No
energy repair, current/density floor or hidden physical resistivity exists.
Massless MHD and fully filled periodic/conducting reference boundaries are
explicit; this is not nonlinear Hall or the free-vacuum formation interface.

All 33 declared cases completed. Exact rotation, pinch and combined equilibria
have maximum spurious radial RMS speed 3.99221e-16 m/s. A finite-amplitude
rotating entropy wave tests actual nonlinear angular-momentum/energy transport;
its nz=192 reconstructed density/rotation relative volume errors are
0.000126545/0.000125193. Torsional-wave reconstructed rotation error is
0.000462097; amplitude controls independently separate quadratic feedback from
spatial error. Stronger contrast and first-order diffusion controls remain
explicit, and no single stronger-contrast run is declared converged.

Independent conducting-Bessel decay reaches nr=128 (6689 accepted steps):
reconstructed toroidal field error 3.35910e-6, first-annulus g error magnitude
1.15482e-6 and physical heat relative error 2.07328e-5. Cell-average error is
separately recorded (1.95079e-6), not confused with field reconstruction error.
Physical quadrature orders 6/10 change reported relative norms by at most
8.21825e-17. Density time-halving ratios are 4.1723/3.9905 and angular momentum
3.7153/4.3801 against a finer same-grid reference, with limiter effects explicit.
Final mass/total-energy relative changes remain within 4.52346e-15/1.13086e-15.

All cases were regenerated under the final fixture/measurement contract, which
honors nondefault gamma and checks continuum fields between cell centers.
Frozen production/measurement SHA256 is
`24ce9b8cdb121167d78ab2ee75792bb9ace362695dd852cbc2384065f4dc69f2`.
Complete source/config/archive/history checks and scientific visual inspection
are recorded in TOROIDAL_MHD_VERIFICATION and `outputs/toroidal_verification`.
All solver, measurement, test and report jobs are terminal at this checkpoint.

Next implement NONLINEAR_HALL_PLAN: compatible nonlinear Hall induction and
Poynting work on the new compressible state, using independently derived
finite-sound-speed cylindrical modes and nonlinear transfer checks. Electron
pressure/thermal and finite-mass current transport, heated-edge compatibility,
spatial coil work, actual vacuum/chamber ends, shock/remapping robustness,
kinetic ordering, converged formation and the experiment UI remain required.
The full objective remains active and the original argon target is unchanged.

## Historical completion audit at the 224-test checkpoint

| Requirement | Authoritative evidence | Status |
|---|---|---|
| Conservative dynamic core | mhd and toroidal modules, independent tests, MHD/positivity and TOROIDAL_MHD_VERIFICATION reports | Bounded nonlinear poloidal/toroidal core implemented; positive traces without cell floors; full shock/update robustness and cold-flow accuracy remain bounded by tests |
| Physical source/plasma/vacuum coupling | Fixed-chamber, boundary, moving/resistive and SURFACE_VERIFICATION reports; SURFACE_COMPATIBILITY audit; Stage 2 guard still raises | Partial: moving resistive free surface and explicit surface refinement are implemented; conditional heated-edge incompatibility, unresolved flow/current, spatial coils, mesh robustness and actual ends remain |
| Regime-appropriate transport / Hall physics | Conditional Coulomb/local regime module, finite-current-inertia component and linear cylindrical Hall/toroidal reference; LINEAR_HALL_VERIFICATION and EXTENDED_PLASMA_PLAN | Partial: linear Hall/finite-mass waves and nonlinear compressible toroidal MHD verified; nonlinear Hall coupling, chosen formation closure, moving electron/ion thermal dynamics and physical edge/kinetic treatment remain incomplete |
| Full independent numerical verification | Prior dynamic/shock/vacuum/surface/current studies; conditional Coulomb checks; 49-case cylindrical Hall/16-mode parallel and 33-case nonlinear toroidal studies with independent species, transport, time and physical-field/energy verification | Partial; 224 passing tests; known false reversal on unresolved current layer, substantial thin-layer B error, moving-flow/cold-pressure/shock robustness and physical surface closure remain |
| Converged FRC formation | Only older stationary-conductor results, unconverged topology duration | Not demonstrated |
| Reproducible experimental interface | Stage 0/1 CLI/playback; dynamic and coupled benchmark scripts, histories and metrics | Partial; no coupled dynamic experiment UI |

Beads issues `frc-goal`, `frc-vacuum`, `frc-transport`, `frc-accuracy`, and
`frc-formation` preserve the remaining scope; `frc-mhd` records the completed
bounded reference core milestone. All numerical jobs at this checkpoint are complete. The optional request for target plasma parameters is not a blocker
to numerical verification or development.

The full objective remains active until all completion evidence is present.

## Nonlinear Hall spatial checkpoint: progress, failed axis gate retained

The full suite passes **244 tests (84.28 s)**. New independently checked
finite-sound-speed Hall modes derive all three branches from compressible
momentum, pressure and Faraday; original electron/ion and Maxwell residuals,
physical energy, conducting/axis traces and limiting dispersions are verified.
The declared reduction has massless cold electrons and adiabatic ion pressure.
A separate finite-amplitude variable-density analytic snapshot verifies local
Hall current, induction and Poynting identities before judging numerical work.

The experimental HallCurl/HallInduction/HallPoynting components cancel assembled
Hall magnetic work through a common electric quadrature and an exact adjoint
of the existing magnetic energy. Both poloidal and toroidal rates are present;
there is no carrier floor, Hall Joule heat, extra fluid force or energy repair.
They are **not enabled in ToroidalMHD**, whose prior verified path is unchanged.

All 42 instantaneous spatial cases completed. Nonlinear nr=8/16/32/64/128
physical magnetic-rate errors decrease 0.259461/0.132568/0.0666923/0.0333967/
0.0167042, and current error reaches 0.00235369. However, first-annulus gdot
error approaches 0.214581 instead of zero. A shrinking physical axis volume
hides that coefficient defect in global norms. This explicitly fails the Hall
integration gate; global conservation is not accepted as physical accuracy.
Local energy-transfer defect decreases to 0.000109583 relative but is not
installed as heating or used to repair magnetic energy.

Maximum assembled Hall work ratio is 2.72060e-17. Independent archive energy
replay backward error is 2.60690e-17 using the actual summation scale. Reference
quadrature 10/14 changes cell comparisons by at most 4.75452e-16; common
quadrature 5/7 changes the nr32 physical induction error by 1.97065e-15.
All three finite-pressure branch/sign Hall-rate errors refine to 3.52–4.00%
at nr32. These are spatial rates, not time-evolved waves. Direct quadratic
Hall amplitude powers are 1.999999987 and 2.000000001; the ill-conditioned
opposite-full-field subtraction and its small-amplitude failure remain saved.

Frozen production/measurement SHA256:
`35cea1122c521983ac6808a8d8f6d27c3cf430428561e2f0de45819598a93996`.
Every declared case/config/source/archive identity passed audit; both final
scientific figures were visually inspected. All test, measurement and report
jobs are terminal. Details: COMPRESSIBLE_HALL_REFERENCE, HALL_WORK_NUMERICS,
HALL_WORK_VERIFICATION and `outputs/hall_work_verification`.

The next action is to fix regular-axis polynomial current/electric/energy
compatibility, then local Poynting transfer and stable nonlinear Hall time
evolution. The original SSPRK2 fluid method cannot simply receive an undamped
imaginary-spectrum Hall RHS. Full electron/ion thermal and finite-mass transport,
heated-edge/vacuum/real-end closure, spatial coils, shock/remapping robustness,
kinetic ordering, converged formation and interactive experiments remain open.
No formation claim or original-argon target change is made. The full goal is
active; this is progress rather than completion or an external blocker.

## Compatible nonlinear Hall checkpoint: verified progress

The full suite passes **257 tests (88.52 s)**. The corrected electric/magnetic
spaces commute under physical Faraday in x=r²/2,z, with explicit compatible
magnetic energy carried through pressure recovery, fluid stages and diagnostics.
The old cell-energy path and failed Hall operator remain reproducible controls.
There is no fitted axis coefficient, density/current floor, energy reset or
extra Hall force. A nonlinear implicit midpoint couples Hall and positive
resistivity to ideal compressible fluid half steps; failed trials are atomic.

All **54 main plus 9 supplemental cases** completed. Axis gdot error falls
0.0203290→0.00526712→0.00132551→0.000331905→0.0000830089 for nr8→128, replacing
the historical nonzero 0.214581 limit. Uniform/strong-density, zero-component,
quadrature and radial-only controls retain convergence. Global Hall work ratio
is at most 4.41337e-18; independent physical energy replay is within 7.94556e-16
of its absolute work scale. Local energy-transfer defect reaches 9.59915e-5
relative and converges near second order, but is not zero or repaired as heat.

The actual coupled solver reproduces all three finite-pressure Hall branches.
At nr32, physical magnetic errors are 3.9960%/4.4048%/3.9631%; initial field
representation already contributes 3.9807%/4.3985%/3.9631%. This is a remaining
first-order spatial field limitation, not accurate-everywhere B inferred from
small budgets. Slow/intermediate/fast velocity errors are
1.0118%/0.2128%/0.2759%. Slow duration is 0.5 radian, the others a full period.
Slow mode retained amplitude improves 0.6961→0.9679→0.9969, so coarse damping
is substantial and explicit. A 100-fold smaller amplitude hardly changes its
nr16 error. Both propagation signs, Hall-free, nonlinear amplitude and solver
tolerance controls are archived; physical and coefficient errors stay separate.

Fixed-grid time-halving magnetic error ratios are 3.94368/3.98028 for the wave
and 3.94105/4.02325 for nonlinear evolution. The finite-amplitude 1-microsecond
flow comparison against restricted nr64 has nr32 magnetic/momentum/density
differences 0.2432%/0.3791%/1.2872% of the fine reference's **change**, rather
than its large background. This does not make nr64 an exact continuum oracle.
Independent conducting Bessel decay has nr64 field/axis/physical-heat errors
1.21612e-5/9.87096e-6/1.33088e-4, with time controls separating spatial error.

All evolved histories conserve total energy/mass to at most
6.63912e-16/5.40031e-16 relative. Independent initial/final energy replay, actual
terminal times, positivity, physical heat monotonicity, source/protocol identity
and lossless archive hashes pass audit. Angular momentum uses an absolute
normalization where its mean is zero. All three scientific figures were visually
inspected. Frozen engine/reference/measurement source is
`9c14e7c757b6d3a51bc0988faff5ff704cd58e2abb8b5ebaada646de7e4d4dd8`;
supplemental protocol and report hashes are separately recorded. All numerical,
report and test jobs are terminal. Details: COMPATIBLE_HALL_NUMERICS,
COMPATIBLE_HALL_VERIFICATION and both `outputs/compatible_hall_*` directories.

### Current completion audit and next action

| Requirement | Current status |
|---|---|
| Nonlinear compressible Hall core | Implemented and independently verified in a bounded conducting periodic reference; first-order physical B representation and limited shock/update scope remain explicit |
| Electron/ion thermal and finite-mass closure | Incomplete: cold massless electrons, single represented gas, prescribed eta; physical heating partition, electron pressure/enthalpy and nonlinear drift-energy terms still required |
| Plasma-vacuum, edge, ends and sources | Earlier moving/free-vacuum references retained; heated-edge incompatibility, spatial coils, real ends and unresolved flow/mesh errors remain |
| Converged FRC formation | Not demonstrated; initialized equilibria and waves are verification fixtures |
| Interactive physical experiments | Existing Stage0/1 interface only; guarded Stage2 and dynamic experiment UI remain required |

Next implement TWO_TEMPERATURE_PLAN's **discrete conservative charge flux and
electron pressure/enthalpy pairing**, then warm electron/ion evolution and
independent collision exchange. The weak work current is not automatically a
pointwise divergence-free cell particle flux; its dual identity must be mapped
and verified before thermal advection. The electron enthalpy drift is required
alongside pressure Ohm/Poynting terms. Finite-mass convective current, drift
energy and momentum stresses must follow from the same species equations.

Return thereafter to original heated-edge/vacuum/real-end physics, spatial coil
work, shock/remapping robustness, kinetic ordering and converged formation.
No easier or colder fixture replaces the original argon target. Beads goal,
transport, vacuum, accuracy and formation issues stay open. The full objective
is **active**; this checkpoint is progress, not completion or an external block.

## Conservative electron pressure/entropy checkpoint: verified progress

The full suite passes **266 tests (88.54 s)**. `DualCurrent` now converts
toroidal magnetic work into exact locally conservative integrated charge fluxes
on full-domain nodal dual volumes. Axis, periodic seam and physical wall trace
are explicit. An independent uniform-current control shows why omitting a
nonzero wall Btheta trace inserts a false return-current sheet and reverses
the last dual-cell current. General electric work includes wall Poynting power.

`ElectronPressure` pairs electron entropy drift by -J/(ne e) with pressure-driven
induction, retaining the enthalpy/pressure energy exchange. Stable logarithmic
spatial means and a magnetic-midpoint/thermal-logarithmic temporal state preserve
both total magnetic/electron energy and discrete physical logarithmic entropy.
The latter is distinct from a convenient linear entropy proxy. Constant-density,
isothermal and adiabatic barotropic fields provide independent zero-battery
limits. Exact gradient identities reduce axis cancellation without a fitted
multiplier, state repair, artificial heat or temperature floor.

All **46 declared records** complete: 17 spatial cases, 28 successful evolutions,
and one deliberately strict nonlinear failure. Pressure-battery gdot errors
decrease 0.0252946→0.00636035→0.00159239→0.000398242→0.0000995696 for nr8→128;
first-annulus error reaches 9.79904e-5 and thermal-rate error 4.10757e-4.
Current face-integral and nondefault gamma/strong-density controls converge.
Uniform-density induction is exactly zero; the isothermal nr64 residual is
1.31973e-9 T/(m s) absolute, with no invented relative zero-reference norm.

Actual electron pressure translates at the analytic electron velocity through
one full periodic transit, with stationary bulk and unchanged Btheta. Physical
continuous pressure errors at nz16/32/64/128 are
0.167636/0.0425302/0.0106771/0.00267227 relative to the perturbation. Nodal and
between-node errors are separate; quadrature uncertainty is at most 6.16e-17.
Opposite propagation and 70% amplitude controls are retained. Fixed-grid
pressure time-halving ratios are 4.04340/4.19891 for drift and
4.00765/4.04658 for finite-field nonlinear evolution; independent DOP853 tests
separately verify second-order temporal consistency.

Starting from Btheta=0, the isolated pressure battery generates 20.0123 microjoules
of toroidal magnetic energy after 1 microsecond at nr64, with the corresponding
electron thermal-energy loss. Compared to restricted nr64, nr32 magnetic and
pressure differences are 0.1122%/0.5289% of reference changes. These are discrete
coefficient/nodal comparisons, not an exact continuum or formation claim.
The strong-density/temperature nr64 case generates 278.219 microjoules and
retains positive electron pressure. Density-scale and solver-tolerance controls
show changes far below the measured mesh errors.

Maximum charge/work/entropy-rate backward errors are
1.38766e-16/1.80985e-16/1.74751e-16. Maximum accepted-history energy drift is
5.70516e-15 relative; physical logarithmic entropy drift is 3.35616e-12 per
particle, including accumulated nonlinear tolerance over long drift runs.
The strict nr64 tolerance 2e-13 control fails at zero accepted time with residual
9.95419e-13; budgets 80/160/320 retain that stagnation and exact untouched inputs.
No failed trial is relabeled as evolved data. The default tolerance is 1e-10.

The full protocol/configuration/source/archive/history and independent physical
pressure/energy/entropy audit passes. Both scientific figures were visually
inspected. Frozen engine/reference/measurement fingerprint:
`5cccf843cfe6ebbf6005d6b87261db7e9ae9decc1277a624b2ec90de8bbc78a7`.
Report-source and figure hashes, supplemental precision controls and lossless
fields are retained in `outputs/electron_pressure_verification`. All numerical,
test and report jobs are terminal. Details: ELECTRON_PRESSURE_NUMERICS and
ELECTRON_PRESSURE_VERIFICATION.

### Remaining integration and completion audit

This new component has **fixed density and no moving ions or Hall J cross B**.
Its constant toroidal-field/density trace boundary generally permits local
pressure Ez/Poynting flux; it is not the old Hall solver's conducting Ez=0 wall.
It is therefore deliberately not inserted as an inconsistent warm closure.
The existing compressible Hall path and all its original limitations remain.

Next implement TWO_TEMPERATURE_PLAN's conservative primal-to-dual mass overlap
and bulk flux map, proving that its divergence commutes with the existing
finite-volume continuity update. Then couple dual electron entropy to changing
density and primal total pressure/ion energy with explicit local enthalpy and
Poynting work, positive species partition and unified boundary assumptions.
Warm coupled waves, collision exchange/resistive heating partition, conduction/
ionization and full finite-mass current transport remain required. No local
energy reset may stand in for that derivation.

The physical heated edge, exact vacuum, real gas ends, spatial coil work,
shock/remapping robustness, kinetic ordering, converged FRC formation and dynamic
experimental UI are still incomplete. The original low-pressure argon target
is unchanged; Stage2 remains guarded. The full goal and its transport/vacuum/
accuracy/formation issues stay active. This checkpoint is concrete progress,
not completion or an external blocker.

## Moving ideal electron/ion coupling checkpoint: verified progress

The full suite passes **283 tests (94.53 s)**. `DualMass` maps actual primal
mass and Rusanov fluxes to the electron dual grid with a locally commuting
divergence, explicit axis/wall half support and periodic seam. Independent
affine physical flux integrals and random open-boundary tests verify this.
The Q1 scalar average is the volume adjoint of the overlap, giving an explicit
electron energy representation. No separate electron density is evolved.

`AdiabaticSpeciesMHD` now evolves moving ions and electron entropy with the same
two ideal RK stages. Changing n supplies electron compression work; positive
ion pressure is recovered from conservative total thermal energy. Equal gamma
keeps total fluid dynamics identical to the existing ideal solver, verified
bit for bit at both orders. No total-energy reset, temperature floor or hidden
retry is used. Resistivity is rejected pending physical heating partition.

All **36 declared evolutions** finish: 30 main and six stronger-gradient/common-
time controls. Independent translating contact and finite-amplitude compressing
simple waves satisfy continuity, momentum, total energy and both species
equations at two gammas. The finest continuous electron-pressure errors are
0.00561259 after one contact transit and 0.000592550 after half a compression
transit at nz256. Finest reconstructed ion-pressure compression error is
0.00150679; corresponding electron/ion temperature errors are
0.000670223/0.00343809. Node and physical-field norms remain distinct.
Pa/eV values and quadrature controls are preserved in the report.

The strong contact has a **negative physical electron entropy change**. The
reported (gamma-1) delta Se/(N kB) decreases in magnitude from 5.83338e-4 to
1.09066e-4 to 2.39562e-5 for nz64/128/256, while physical pe errors decrease
0.152234→0.0417378→0.0123756. This scheme has no proven entropy inequality or
shock heating partition; exact total energy does not remove that limitation.
Strong Mach .12 compression pe errors reach 0.000738578. Positive accepted and
reconstructed species states are retained throughout all recorded runs.

Fixed-grid time ratios are 4.03370/4.50765 contact and 4.10340/4.30592 compression.
Half-CFL changes are 0.401%/4.19% of the base node errors. A radial refinement
control's apparent 5.14219e-6 pressure difference becomes 1.61383e-15 when time
steps are matched, identifying the original CFL confound. Reference-density
scale changes pressure by 1.06722e-15. Opposite translation, gamma 1.4, first-order
and nearly cold-electron controls finish. Physical-norm quadrature uncertainty
is at most 2.41243e-10 absolute in relative-L2 units.

Maximum accepted-history mass/total-energy/linear-proxy drifts are
5.40031e-16/8.87805e-16/1.12490e-15 relative. Independent physical magnetic energy,
quarter-volume species density, electron/ion energy and logarithmic entropy
replay pass. Frozen engine/reference/measurement fingerprint:
`0e37c39afe2dc22bd122cd6d21ec5e21981cfa97b7d5d59e154bd263d88c73b4`.
Main and supplemental source/configuration/archive/history audits pass; both
scientific figures were visually inspected. All jobs are terminal. See
WARM_ADIABATIC_NUMERICS and WARM_ADIABATIC_VERIFICATION.
All 36 cases were regenerated under this final source after compacting history
storage into lossless columns. Every accepted value, endpoint, physical metric
and NPZ hash matched the earlier layout exactly; both figure hashes are
unchanged. Saved collections shrink from about 143 MB to 43 MB without losing
diagnostics. The earlier c4dca9... working fingerprint is superseded.

### Remaining integration and completion audit

This is a **common-velocity ideal thermal reference**, not a coupled warm Hall
solver. Its ideal electric wall trace does not validate the full pressure
electric field along a conducting wall, even in curl-free one-dimensional
fixtures. The isolated pressure/current and cold Hall components remain separate.
The earlier 63 Hall-case source hashes were checked against 16969a9 after the
new mass-flux interface changed source identity; they were not relabeled as
new-source records. The isolated electron fingerprint is unchanged.

Next derive and implement local conservative electron-current enthalpy and
Poynting transport on this moving state with mutually consistent magnetic
boundary traces, then verify warm Hall waves and nonlinear evolution. One
candidate worth deriving is natural fixed Btheta(R)=0 with free pressure Ez(R),
which makes boundary magnetic power zero without requiring constant wall n;
it is not the previous conducting Ez=0 boundary. This candidate is not yet
implemented or verified. No mapped local energy reset may replace the flux.

Then physical collisions/heating partition, conduction/ionization and full
finite-mass current/drift energy must follow. Original argon validity, heated
edge, exact vacuum, real gas ends, spatial coils, shock/remapping/entropy
robustness, converged FRC formation and dynamic experimental UI remain required.
Stage2 stays guarded; the full goal and transport/vacuum/accuracy/formation
issues remain active. This is progress, not completion or an external blocker.

## Warm Hall coupling and pressure-split limitation checkpoint

The full production suite passes **293 tests (94.10 s)**. A zero-Btheta/free-Ez
magnetic trace now supports full pressure gradients and varying wall density,
with exact curl-gradient/dual-current work identities. The moving warm reference
adds conservative pressure/Hall Poynting power and electron-current enthalpy
face fluxes. No local total-energy reset, pressure floor or hidden retry is used.
The five-stage coupling is experimental and retains eta=0.

All **58 main cases** finish and are audited. Conservative local pressure-energy
rate errors decrease .0940091→.000395013 over nr8→128; local ion-work residuals
also decrease. Warm branch2 electron/ion pressure errors at nr32 are
.00347082/.0173594, while physical magnetic error remains .0394869. The slow
branch is still inaccurate: finest electron/ion pressure errors are
.164298/.829787 and ion-temperature error is 2.34858 of its small perturbation.
Half-CFL and tenfold amplitude reduction do not remove the slow error.
Fixed-grid pe time-halving ratios are about four for isolated pressure, fast
warm waves and nonlinear coupled evolution. Nonlinear nr16 versus nr32 physical
pe/B differences are 4.19%/9.99% of the reference trajectory changes, not a
converged formation claim. Maximum main mass/total-energy/proxy drifts are
5.41e-16/2.15e-15/5.82e-15. Independent physical energy, density, species and
field-norm replay pass; all three scientific figures were visually inspected.

The boundary resistive mode retains a **-15.21% instantaneous first-annulus
g-rate error**, even though its evolved coefficient errors decrease. Evolved
time controls show finite time error too. This is a separate unresolved axis
consistency issue; resistivity is not enabled in the warm model.

Five supplemental runs finish, but the sixth exposes a stronger limitation.
At nr64 the maximum local pressure-only ion-energy error grows to .0921 Pa
despite a decreasing volume RMS; finer time steps retain about .0914 Pa.
At nr128, three dedicated replays with 320/640/1280 planned steps reject
negative ion pressure at about .95 microseconds, at cell [0,100]. Electron
pressure stays above 51 Pa, nonlinear residuals are below 1e-10 and global
energy drift is at most 2.08e-16. Accepted and rejected states are separately
archived; none of these failed trajectories is labeled completed.

An independent local continuum pressure-only symbol has discriminant
-7.40906e7 m²/s² at an interior point of the original nonbarotropic reference,
including background-current advection. Its frozen short-wave growth is
4303.79 |k| per second when |k| is in inverse meters. This rules out assuming
general hyperbolicity of that isolated split. It is not a proof that full warm
Hall is illposed or that this interior symbol caused the observed axis failure.
See WARM_PRESSURE_SPLIT_LIMITATION for the derivation and preserved evidence.

Frozen engine/reference/measurement fingerprint:
`ff529e0fab60fb0888806589056142c7040b3c66fcff0cfde669635e1fbb7c20`.
Complete data consist of 58 main completions, five supplemental completions and
three rejected replays, plus the analytic symbol scan. Main/supplemental source,
configuration, archive, accepted-history and failure audits pass. Reports hash
their helper sources and inspected figures. Earlier warm-ideal source is
unchanged; old electron records were source-verified at419c951 and old Hall
records at16969a9. They remain historical measurements of those sources.
All numerical/test/report jobs are terminal at this checkpoint.

### Immediate continuation and remaining objective

First analyze the combined nonideal principal behavior, then derive and verify
a joint Hall-plus-pressure stage with one magnetic midpoint and thermal mean,
shared current work, conservative face energy and atomic validation. Joint
fixed-density splitting is a candidate, not a proven cure: degeneracy for
k perpendicular to B may require retaining ion/density or finite-mass dynamics.
Preserve the exact failed fixture and time/axis/maximum-norm controls; do not
add a convenient bias, arbitrary heat or numerical atmosphere to hide failure.
Address the slow-wave and standalone resistive-axis errors explicitly.

Then physical collision/heating partition, conduction/ionization and full
finite-mass current transport remain. The original low-pressure argon target,
heated edge/true vacuum, actual gas ends, spatial coil work, shock/remapping/
entropy robustness, kinetic validity, converged formation and interactive
experiments are unfinished. Stage2 stays guarded. The full goal and its
transport/accuracy/vacuum/formation issues remain active. This checkpoint is
progress with a diagnosed limitation, not completion or an external blocker.

## Finite-mass current and precision checkpoint

The separate `CurrentThermal` subsystem now evolves magnetic/current energy
and both species thermal entropies jointly at fixed positive density and
constrained zero barycentric flow. Physical reduced electron/ion mass and full
nonlinear convective inertia are retained; no effective mass, heating repair,
atmosphere or state floor is added. CURRENT_THERMAL_NUMERICS gives parent species,
constraint reaction, boundaries and the full physical energy flux.

All 41 baseline cases are audited: 40 finish and the massless nr128 original
fixture rejects at its 0.75 microsecond trial, with positive ion pressure at
the last accepted state. Physical mass reduces the unchanged original fixture's
local pressure-defect RMS from 0.0591249 to 0.000254304 Pa across nr8→128;
maximum/axis errors also decrease. A finer-time control confirms this particular
error is spatially dominated. Full component/conservative energy drift stays
below 5.80e-16/3.87e-16. Exact nonlinear waves and full face powers refine, but
physical B/J interpolation errors remain roughly first order. A local massless
joint/moving symbol retains nonhyperbolic directions; no global theorem or
claim that it caused the numerical axis failure is made.

The independent intermediate convective-inertia forcing has a nonconvergent
first-annulus error near 0.2354. The separately verified physical inertia
inversion gives full-rate L2/axis errors 1.00341e-4/6.56197e-5 at nr128;
isolated convective L2/axis errors are 2.74914e-4/6.89978e-3. The positive
magnetic-space elimination is equivalent to the original static mixed system,
verified against fresh direct solves and a separate Fourier/Chebyshev continuum
reference. The interrupted direct archives remain labeled partial.

Four uniform-bias failures were reproduced exactly from the frozen source.
Their subtraction of nearly equal weak currents is now replaced by direct
Faraday-increment current evaluation. Small-signal Ohm and stored-state backward
errors are checked separately at unchanged tolerances. Eighteen ten-step
stationary controls and 15 nonlinear/weak-wave evolutions complete and are
independently audited. Maximum stationary B error/bias is 3.54e-14; tiny waves
retain the same physical refinement as ordinary waves. Energy/rate replay
roundoff limits under strong backgrounds are explicit, not hidden by a pressure
repair. The six repeated baseline metrics agree near roundoff.

Baseline source `6c3e50a2d595ee1e2bbe14c5a6679d8da0855ab99143fa802f3422543dec8cd8`
and exact report helpers are retained in a hashed source snapshot. Current
engine source is `80fd45bddda51871f7f36682e8d0bdc4495d0c41430191c25ef03e06ed47c4d2`.
Reports verify every relevant source/archive hash. Baseline convergence/limits
and the physical-rate plot were visually inspected. See
CURRENT_THERMAL_VERIFICATION for the detailed errors and reproduction paths.
The final full suite passes 316 tests in 103.41 s. All local numerical/audit/test
jobs are terminal; the redundant remote direct solve was explicitly terminated
and its process exit verified. Issue state is exported locally; no publication
or remote Git configuration was performed.

Next derive a genuinely moving finite-mass update with drift stress, changing
1/n current kinetic energy, bulk/current transport and velocity-gradient terms
in one conservative contract. First construct an independent manufactured
moving, varying-density reference that checks both parent species and local
energy. The fixed-n component must not be bolted into the moving split by only
changing its density matrix. Physical collision/heating partition, transport/
ionization, the original heated free edge and exact vacuum, real ends and
spatial coil work, shocks/entropy/kinetic validity, converged formation and
interactive experiments remain open. Stage2 remains guarded, the overarching
FRC goal remains active, and there is no external blocker or repository remote.

## Moving finite-mass rate and continuum energy checkpoint

The full suite passes **322 tests in 107.36 s**. `MovingCurrent` now includes
bulk B* advection, velocity curl, actual prescribed n_t/n and kinetic/mixed
scalar potentials with physical reduced mass. A positive magnetic-space solve
handles both poloidal and toroidal instantaneous rates. No density, momentum
or thermal trajectory is advanced by this component.

Independent manufactured parent species check continuity, momentum, Ohm,
Faraday, the B* identity and direct full energy flux at physical Ar+ masses,
the original n0=1e20 regime, equal masses and consistent massless limits.
Pointwise identities agree within 5.09e-15. Local energy differentiation
refines approximately fourth order; dropping mixed bulk/current transport
retains a finite 9.68227e-4 relative defect in the n0=1e18 argon fixture.

All **19 rate cases** finish and pass source/configuration/archive/physical
audits. For n0=1e18, nr8→128 B_t errors decrease 0.217882→0.0116252,
J_t 0.275820→0.00459108, E 1.54929→0.00643579 in physical volume L2.
Finest first-annulus errors are 0.000283210/0.00711149/0.00708714.
Volume B_t and axis J_t/E are approximately first order; axis B_t is second
order. Original-density and phase/mass/quadrature/tolerance controls complete.
The massless nr64 J_t volume error remains **0.403719**, with latest observed
order only 0.588; this numerical limitation is retained for investigation.

Removing density/vorticity/bulk terms separately changes the nr32 magnetic
rate by 1.89967%/0.146032%/28.0284%. Maximum full-E Faraday replay error is
1.3718e-11; the changing-density physical energy-derivative balance is at most
1.3690e-10 of component powers. These are instantaneous identities, not a
conservative complete-update energy theorem. Higher norm quadrature changes
relative L2/axis metrics by at most 1.33e-15 absolute.

The initial nr128 residual rejection and target nr32 energy-cancellation audit
failure are preserved with their exact source. Explicit linear refinement of
the original residual and stable polarization of the unchanged energy functional
resolve them without relaxing tolerances or repairing physical states. One
linear refinement is used at nr128 and in the tighter nr64 control. All final
cases share fingerprint
`d3984fb056edb1aa150d030937b9bce6c21d292b84ef80305b19b52f88af5489`;
the previous current/thermal fingerprint 80fd45b... remains unchanged.
The final source snapshot, report/helper/field hashes, inspected figure, full
test log and failed-attempt audit are retained. All jobs are terminal; no
remote Git, publication or background numerical job was introduced.

Next derive and implement a **discrete conservative moving finite-mass update**
with rho=(mi+me)n, actual continuity flux/n_t, a consistent regular velocity
map, drift-stress momentum/angular momentum, species thermal transport and full
bulk/current/Poynting/enthalpy energy flux. Verify local pressure/energy and
independent finite-mass compressible waves/manufactured moving trajectories
with spatial/time convergence and atomic failure. Do not bolt the rate onto
the old split. Physical collisions/transport/ionization, heated free edge and
exact vacuum, real ends/spatial coils, shocks/entropy/kinetic validity, converged
formation and interactive experiments remain open. Stage2 stays guarded and
the overarching goal remains active. See MOVING_CURRENT_NUMERICS/VERIFICATION.

## Conservative moving finite-mass integration checkpoint

The full suite passes **328 tests in 120.97 s**. FiniteMassMHD now evolves
rho=(mi+me)n, all momentum components, conservative total energy, psi/g and
electron entropy through four complete unsplit RK4 derivatives. The actual
mass flux supplies dual n_t, induction retains physical reduced mass and bulk
velocity gradients, and shared fluxes carry drift stress, species enthalpy,
mixed current energy and full computed Poynting power. All stages validate
positive species states before atomic commit. No physical step retry, thermal
floor, artificial mass, atmosphere or total-energy reset is introduced.

All **21 complete-rate cases and 18 declared trajectories** finish and pass
source/configuration/archive/physical audits. Parent manufactured species and
independent exact-mass compressible waves verify rates and final states. The
manufactured and unforced nonlinear runs last 50 ns; the fast waves last one
quarter period (29.1294 ns). Slow/medium branches have rate evidence only.
Global mass/energy-budget errors are at most 4.44e-16/3.26e-16; independent
physical energy replay differs by at most 2.55e-14. Explicit source work matches
independent parent total-energy change within 5.08e-12 J. All accepted species
pressures are positive. These conservation properties do not prove local accuracy.

The finest manufactured nr128 mass/total-energy rate errors are
.000861842/.000529452, but pe/pi rate errors remain .00774879/.143364.
The pi-rate maximum increases slightly from 30535.9 to 31209.6 Pa/s at nr64→128,
despite decreasing volume/axis errors. The slow-mode nr64 pi-rate error is
41.4570 times its physical signal. Tenfold amplitude reduction leaves these
errors essentially unchanged. Preserve maximum norms as well as L2 evidence.

At nr32, evolved manufactured ion pressure has .00362760 relative error against
the full pressure and **.783580 against the true pressure change**. The original
n0=1e20 case retains .00250372/.540816. Fast-wave magnetic/electron/ion pressure
perturbation errors are .0398052/.0105524/.595516. The weak-wave control preserves
that limitation. Fixed-grid dt halving has a much smaller effect than mesh error.
Magnetic time ratios are about 15.4; pressure ratios range 3.15–5.19, so formal
RK4 order is not asserted for every quantity. Unforced electron physical log
entropy decreases by 5.13e-6→3.95e-7 per particle under refinement; the method
has no entropy inequality or physical shock partition.

The slow-wave continuity error is chiefly the limited central reconstruction,
not Rusanov density diffusion. At nr32 those error/contribution norms are
28.6816 and 1.03981 times the true rate. A separately tested cubic reconstruction
of exact annular linear momentum averages reduces the density-rate error to
.0202495 at nr32 and .000179355 at nr128. With actual finite-amplitude momentum
it instead levels near .00311608 against the linear oracle. This candidate is
**diagnostic only**, not adopted into the coupled solver. Correcting mass flux
alone would leave enthalpy/kinetic/momentum work inconsistent.

Baseline rate source is 8937eefb61c197962881e7d844bee6d9f5dfce45aff93025c1c35e00dba63ea7;
expanded trajectory source is 3edae52928f5cdaa20ca240f11d6fa729a137be327c9b1ffc64f09c6440f42fb.
Baseline production sources are unchanged; only provenance coverage expands.
Prior moving-rate d3984fb... and current/thermal 80fd45b... fingerprints remain
unchanged. Full records, fields, histories, physical/source audits, diagnostic
arrays, source snapshots, inspected accuracy figure and test log are retained
in outputs/finite_mass_initial_attempt and outputs/finite_mass_evolution.
See FINITE_MASS_NUMERICS/VERIFICATION for scope and complete measured errors.

Immediate next work is to derive shared higher-order particle, momentum,
enthalpy and kinetic-energy transport. Verify a translating two-temperature
contact and all finite-mass wave pressure rates, including near-solenoidal
cancellation, before adopting a reconstruction. Retain independent parent
pressure/energy, axis/max norms, exact uniform flow and atomic failure checks.
Then revisit physical collision/heating partition, conduction/ionization,
original heated free edge and exact vacuum, actual gas ends/spatial coils,
shocks/entropy/kinetic validity, converged formation and the experimental UI.
frc-7mm and the overarching goal remain open; this is progress with a diagnosed
accuracy limitation, not an external blocker or a completed formation simulator.

## Shared smooth transport and magnetic cancellation checkpoint

The full suite passes **336 tests in 118.11 s**. The separate SmoothFiniteMassMHD
reference uses shared cubic conservative face transport, mean-preserving
transverse quadrature and a consistent ratio of reconstructed dual U and n for
entropy flux. The historical limited solver is unchanged. This central smooth
method introduces no physical/artificial diffusivity, floors, hidden fallback
or atmosphere. It has no shock/entropy/positivity theorem; invalid stages reject
atomically. Exact vacuum is still an unfinished coupled-model requirement.

All **31 rates and 27 trajectories** pass the independent source/field/history,
parent-reference, physical component-energy/entropy and conservation audit.
Energy-budget residual is at most 3.36e-16; normalized integral drift is at most
4.97e-15. All measurement and audit processes are confirmed terminal. Source
dd060015dc510d312f22a53588d2809ad806de02de1dadb4db956e8887ab4ead and its snapshot
are in outputs/smooth_finite_mass_verification. Earlier finite-mass 3edae529...,
moving-current d3984fb... and current/thermal 80fd45b... identities are unchanged.
The first smooth attempt, original electric-kernel failure, serialization error
and orphan fields remain preserved. All 49 old complete records and both orphan
arrays reproduce exactly; new contact records supply the missing histories.

At nr128, manufactured pi-rate error improves from .143364 to .0268603; its
maximum now decreases to 3225.44 Pa/s. Slow nr32/64 pi rates improve from
117.248/41.4570 to .0691434/.00883534. Full fast-wave nr32 pi error improves
from .595516 to .0538315; magnetic error remains .0398052. Manufactured nr32
pi error in its small physical change remains .702663. Exact global energy
and accurate thermal rates do not establish accurate magnetic dynamics.

Three unmagnetized two-temperature contacts finish a full 100-microsecond
transit. Finest nz64 nodal pe/cell-mean pi errors are 6.63e-6/.000316986;
continuous pe/pi errors are .0000750472/.0110717 and temperature errors
.000498286/.0127339. Initial projection/representation errors, reverse flow,
limited entropy and half-dt controls are explicit. Small false B/J drift is
not monotonically convergent. Three magnetized contacts last only 50 ns and
must not be described as full-transit magnetic verification. Physical entropy
defects remain and reverse sign with flow in the contact control.

The exact separable electric-kernel representation removes a contact roundoff
rejection while preserving strict small-signal residuals and original full
equation backward errors. The full recombined Faraday error relative to its
tiny rate remains explicitly .016679 in the preserved example; there is no
claim of well-conditioned full-array reconstruction or arbitrary-kernel support.

Slow magnetic-rate errors remain 1004.98/263.069/66.8838/16.8431 at nr8→64;
slow conserved-energy-rate errors remain 106.207/26.5782 at nr32/64. Exact
velocity/current forcing controls identify incompatible lower-order field
representations, amplified by roughly 19200-fold ideal/Hall cancellation.
Neither mesh refinement alone nor fitted velocity from numerical current is an
adequate resolution. Next derive shared higher-order velocity and magnetic/
electric bases with physical Hodge/inertia quadrature and commuting curl.
First verify independent polynomial fields, Ampere/current convergence, all
wave branches and full original Ohm/Faraday/energy contracts, then repeat the
complete smooth trajectories without changing their physical inputs.

The numerical/physical requirements in frc-7mm, frc-accuracy, frc-transport,
frc-vacuum and the overarching goal remain open. Continue collision/heating
partition, conduction/ionization, heated free-edge exact vacuum, real ends/
spatial coils, shocks and kinetic validity, converged formation and the UI.
No FRC stability or particle-confinement claim follows from this checkpoint.

## Higher-order compatible field and harmonic precision checkpoint

The full suite passes **343 tests in 118.34 s**. New `compatible_polynomials`
and `compatible_geometry` supply local degree-1/3/5 nodal/interval-average bases
with commuting derivative, unchanged Faraday incidence and physical tensor
electric/magnetic masses. `compatible_velocity` reconstructs actual annular
conservative averages. `compatible_current` uses physical finite species masses
and constant-background inertia. Five independent polynomial/integral/current/
energy/Faraday/mixed-system tests verify the component, not just its algebra.

The first 25-case rounded-total study exposed a precision limit: slow nr32/64
B_t errors .0117924/.0245827 despite far smaller spatial error in the weak
perturbation-only control. A hundredfold weaker full-bias wave gave 1.17396.
This complete attempt is preserved under source 485edaa1fdb9bffcb76b2a4e9196e3e4b1ddc071b20f319b4c3062b42ebae667.
An independent three-way probe distinguishes current-operator cancellation from
digits already lost in storage; subtracting bias afterwards does not recover them.

`HarmonicMagnetic` now stores an explicit constant physical bias and independently
retained increment. Uniform axial bias has exactly zero current under the stated
fixed-flux periodic boundary. Full field/local energy keep the bias and local
cross term; its cross integral is zero by exact boundary values. Two polynomial
tests verify weak currents through amplitude 1e-12 beside that unchanged bias,
small-increment energy and independent full local energy. No signal threshold,
fitted velocity/current/mass, energy repair or numerical atmosphere is introduced.

All **25 affine cases and 25 rounded comparisons** pass the independent full
source/input/field/physical-norm/original-equation audit. Main full-bias degree5
slow B_t errors at nr8/16/32/64 are .789352/.0141969/.000190057/2.36657e-6.
Finest intermediate/fast errors are 6.68e-9/2.51e-9. Slow current-rate error
decreases to 7.73307e-5, less favorable than magnetic accuracy; axis/max errors
are retained. Hundredfold weaker signal, opposite propagation, original n1e20
and quadrature controls retain the improvement. The exactly matched earlier
linear diagnosis used the same physical inputs, without analytic forcing fields.

Largest independent physical energy discrepancy is 2.75e-14, including error
normalized by the small increment energy. Original backward errors are below
1.83e-12, but full-E Faraday replay relative to its tiny rate reaches 5.88e-10
after restoring pressure potential. Rounded export still loses current accuracy;
its worst current-relative Ampere replay error is .000129808. These distinct
measures are not mislabeled as small-signal solver residuals. All jobs are
confirmed terminal, and the scientific figure was visually inspected.

Final source is 8cfb1a2f89f7ee58fadae2e47ae9e4fbb5c49037face19da117518188ed8d720
in outputs/compatible_wave_verification; report, snapshots, failures, all arrays,
test/worker logs and checkpoint hashes are retained. The smooth nonlinear source
dd060015... and prior finite-mass 3edae529... remain unchanged. The new component
is **not yet adopted in nonlinear moving evolution**, and no trajectory energy,
entropy, shock, formation or confinement claim follows from these rates.

NEXT derive positive variable-density physical current inertia and local energy
derivatives on this field basis, with actual conservative mass-flux n_t. Derive
the common species transport/pressure-gradient adjoint and full velocity/vorticity,
stress/Poynting/current-enthalpy work before replacing the old state partition.
Retain local bias cross energy and exact physical masses. Verify full independent
parent manufactured equations and coupled wave/contact trajectories, including
small-pressure-change, axis/max, time, quadrature and atomic failure evidence.
The original heated free edge/exact vacuum, transport/ionization, real ends/coils,
shock/entropy and kinetic validity, converged formation and UI remain open.


## Variable-density physical inertia and the axis thermal limitation

This continuation made verified component progress. The final suite passes
**351 tests in123.17 s**, including eight new focused checks. All14 prescribed-
continuity manufactured cases completed and passed original-equation, source,
raw-array and independent physical local-energy audits. Every worker is terminal.

`CompatibleDensity` uses actual conservative cell-average number density and the
same linear reconstruction for supplied continuity n_t. Its physical adjoint
pairs the changing-current-energy density work with the actual cell state.
`VariableCompatibleCurrent` supplies positive variable-density inertia through
tensor products, with physical masses and the unchanged original induction
solve. Local J dot J_t, −J² n_t/(2n), and full B dot B_t (including bias cross
work) remain separate. Negative/nonfinite reconstructions or unrepresentable
inertia reject, without floors or energy repair.

The14-case degree3/5, resolution, quadrature, physical-density and flow-sign
protocol is instantaneous and explicitly manufactured. Finest degree5 B_t/J_t
relative errors are4.96985e-10/4.65041e-10. Local current-energy-rate and its
density-contribution errors are8.93213e-10/1.17057e-9. Independent root-product
field and primitive-polynomial density quadrature agree with production local
energies/rates within4.57e-15. Original backward error is below9.94e-13; checked
solver residual below1.33e-11 at1e-10 tolerance. Full-array Faraday rate-relative
error reaches5.77e-10 and remains separately exposed.

`CompatibleCharge` replaces the old diagonal-g stream with (H_g g)/(dx dz),
restoring graph divergence and pressure/current work pairing. The explicit
`CompatibleSpeciesPressure` bridge retains the old lumped nodal thermal state.
Its species work and fixed-n logarithmic entropy cancel near roundoff, but its
physical pressure E remains roughly second order and entropy rate is inaccurate.
Finest electron/ion whole-domain nodal entropy-rate errors are6.35129%/6.28264%.
The first four radial rows contain99.8072% of the electron squared error; their
local relative error grows toward39.09% with refinement rather than converging.
The interior relative error is0.09767% at nr64. This local axis defect is the
next thermal target; global conservation or shrinking error volume cannot
justify adoption. The scientific figure and radial profiles were inspected.

Final source fc106dc19b9bc208aaf7d9b36a43dc69ea49f4d49d126b4f4dfd2ecd85c59e97
is saved in outputs/variable_current_verification. The complete initial source
310809e2 and all14 cases/audit/test logs remain in
outputs/variable_current_before_range_guard. Final finite-range guards reproduce
every saved numerical array bit for bit across all14 cases. The old smooth
source dd060015 and linear compatible source8cfb1a2f remain unchanged.

NEXT derive a common physical higher-order thermal density/energy/transport
metric and its pressure-work adjoint, explicitly correcting the nonconverging
near-axis lumped entropy rate. Then derive actual conservative velocity/vorticity,
stress and complete bulk/current/enthalpy/Poynting energy transport, and repeat
independent parent-species and full wave/contact trajectories. This component
has not replaced SmoothFiniteMassMHD. Heated free edge/exact vacuum, physical
transport/ionization, real end/coil coupling, shocks/entropy, kinetic validity,
converged formation, three-dimensional stability/confinement and UI remain open.
No remote exists; work and issue state are saved locally without publication.

## Physical thermal metric and full-transit current transport checkpoint

This continuation made verified component progress. The complete suite passes
**361 tests in 121.42 s**. All 18 rate cases and nine prescribed-current
transport cases completed, with original-equation, source/raw-array and independent
physical audits. Every worker, audit and report process is terminal. The
six-panel scientific figure was visually inspected.

`EntropyScalar` supplies the physical number-weighted scalar metric on the same
reconstructed conservative density as variable-current energy. The new nodal
specific log entropy in `CompatibleLogThermal` gives positive representable
pressure, a pressure field paired with physical thermal work and explicit
density/entropy energy derivatives. Current transport conserves physical log
entropy; its exponential proxy is a separately measured approximate invariant.
No pressure floor, current/mass fitting or energy repair is introduced.

Degree-5 electron first-four-row log-rate errors are 1.02400e-3, 2.96189e-5,
5.42997e-7 and 8.84960e-9 at nr=8/16/32/64. The matched old bridge approaches
0.390908 instead of converging. Finest continuous physical proxy-rate error is
1.30658e-8 versus the old representation's 0.117630 in that same norm. Finest
electron/ion pressure-E errors are 4.37e-8/4.21e-8. Independent local energy
agrees within 6.83e-15; pressure work and declared bulk-energy exchange defects
are at most 1.85e-14/1.37e-14. These conservation measures are not accuracy claims.

Constant entropy has exactly zero stored current-driven entropy/magnetic rate.
Uniform temperature at varying density is not exactly representable: its false
B_t RMS converges 0.0416574, 0.000921031, 1.58422e-5, 2.53764e-7 T/s. Restoring
the full rounded E in the finest case gives a Faraday replay difference of
0.00242218 relative in the physical magnetic norm, or 6.15e-10 T/s RMS. Its raw
coefficient-space relative difference is 0.0203616. These distinct conditioning
limits remain explicit despite small original-equation backward errors.

All nine transport fixtures last one complete axis transit, 1.20801 microseconds,
with prescribed density/current and a sheared analytic entropy wave. Full-domain
entropy-signal error is 0.0846755, 0.00203068, 9.31941e-5 at nr=8/16/32; physical
first-annulus error decreases to 9.65e-9. Half-dt changes the nr16 final field
by only 1.73e-8 of the signal. Physical entropy drift is at most 2.22e-16;
proxy/internal-energy drift and initial projection error remain separate. The
finest signed proxy/internal drifts are 3.61e-12/-1.47e-12. This is a prescribed-
current transport verification, not a force-balanced or full moving plasma.

Source b0dba42d2dd23742246776328413c3914ecfebaadf158e72d9b157c10d6485e7 and
complete records are saved in outputs/log_thermal_verification and
outputs/log_entropy_evolution. Original variable-current fc106dc1, compatible
linear 8cfb1a2f and smooth nonlinear dd060015 identities remain unchanged.
See LOG_THERMAL_NUMERICS/VERIFICATION and the hashed checkpoint for full limits.

NEXT derive the common moving state and total-energy contract before adoption:
shared conservative continuity and bulk entropy advection, true conservative
velocity/vorticity, momentum/drift stress, density/pressure work and complete
current/enthalpy/Poynting transport. Resolve the two-nodal-entropy versus old
cell-energy partition explicitly; then verify independent parent-species rates
and full nonlinear wave/contact trajectories. Keep frc-7mm, frc-accuracy and the
overall goal open. Heated free edge/exact vacuum, physical transport/ionization,
real ends/coils, shocks, kinetic validity, converged formation and UI remain.
Work and issue state are local; there is no configured remote or publication.
