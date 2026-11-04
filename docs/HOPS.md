# HOPS — Harmonically Oscillating Plasma Sphere

The user's positive-cage experiment with dynamically timed driving. This first
implementation is one feedback law, not a claim that its thresholds or initial
parameters are optimal. The plasma shape and field evolve freely in full 3D.

## Physical question

Can timing feedback sustain coherent compression more effectively than fixed
pulsing, at comparable measured source work and particle retention?
POPS resonance was observed in an electron-created virtual cathode, not proof
of fusion gain or of an ion well inside an empty positive cage.
[Park et al., PRL 95, 015003 (2005)](https://doi.org/10.1103/PhysRevLett.95.015003).
Electron supply and space-charge neutralization during compression remain
central questions. [Nebel et al., APS DPP 2005](https://meetings-archive.aps.org/dpp/2005/lp1/83/).

An independent Gauss-law derivation for uniform net negative charge gives
omega_i²=Z e²(ne-Z ni)/(3 epsilon0 mi). At fixed density, single-ion harmonic
frequency is independent of amplitude. For a quadratic well,
omega_i²=2 Z e Delta_phi/(mi R_well²); mean-square cloud radius can breathe at
2 omega_i. Fixed charge compressed uniformly would give n proportional to
R^-3 and frequency proportional to R^-3/2. Those assumptions need not hold
with real electron/ion losses and a changing well. No such scaling is imposed
in the actual force or this controller.

## First declared experiment

- Full open 15 cm cage, 2 mm wires, assumed spherical 40 cm grounded chamber;
  absorbing steel and actual-mass 3D Vlasov–Poisson, as in [IEC_3D.md](IEC_3D.md).
- Actual 0–5 µs, 5 ps steps. Positive ceiling +3 kV, 20 ns voltage response,
  common 200 ns on warmup. This is an ideal voltage source, not circuit design.
- After warmup, save the measured live-ion RMS radius R0 about its measured
  centroid. Command off on inward crossing of 0.8 R0, on on outward crossing
  of R0; 10 ns radius sensor, 20 ns minimum dwell. These are controller choices,
  not particle constraints. Sensor filtering does not change fields or charge.
- Actual command onsets define repetition frequency. No contraction or
  re-expansion means no invented successful pulse cycle.
- DC and fixed 1 µs / 50% duty controls follow sequentially, with the same seed,
  source, ceiling, response and warmup. The fixed period is an input, not a
  found resonance. Equal peak voltage does not mean equal delivered energy.
- Seed n=1e11 m^-3, radius 5 cm, Te=1 eV, Ti=0.05 eV, 600 macros/species,
  singly charged mass-2mp ions. Ideal wall electron source 1 mA at 2 eV.
  These are exploratory inputs, not the historical air discharge or an optimum.

Frozen campaign: `outputs/hops_v1/`; actual run:
`outputs/lab_runs/hops_v1_hops/`. The [local viewer](http://localhost:8792/?run=hops_v1_hops)
shows saved 3D positions, radius/voltage, potential and raw budgets.

## Evidence and limits

Two actual source-off +3 kV / 20 ns startups at 10/5 ps completed. Maximum
energy residual divided by particle kinetic/absorbed/injected throughput fell
0.3780%→0.1128%; charge error stayed below 1.87e-27 C. Both final states are
finite. This is timestep evidence; existing FEM small-well force and phase
errors, particle/mesh/evolving Debye resolution and physical losses remain open.

Record centroid/principal radii, live-ion RMS/quantiles, losses, translation,
breathing and residual kinetic energy, actual intervals and full states.
Residual motion includes shear/rotation and is not a temperature. Selective
loss can shrink radius; radius alone does not prove compression. Source work,
electrode/particle charge, field/particle energy and absorption stay explicit.
No imposed well, harmonic electron background, symmetry filter, mass scaling,
charge smoothing or energy repair. Collisions, ionization and fusion are absent.

Current priority (user, Oct 25): measure a pulse and its relaxation before
selecting adaptive timing. The v1 radius threshold never triggered: a controller
that requires contraction before releasing cannot measure relaxation when that
contraction does not occur. A prescribed electrode pulse is an experimental
boundary input; it does not prescribe the plasma's shape or response.

The new single-pulse experiment commands +3 kV for 0.2 µs, then zero through
actual 5 µs, with 20 ns rise/fall response and 5 ps particle steps. The ideal
1 mA/2 eV wall emitter also stops at 0.2 µs, so the later cloud receives no
ongoing particle-source power. All other mesh/seed/species inputs match v1.
Existing fixed-drive support uses a 10 µs period and 2% duty, zero warmup:
no second pulse occurs inside this 5 µs run. This is a single-pulse response
probe, not the claimed adaptive HOPS waveform or an optimized pulse width.
Its config, source and short falling-edge timestep comparisons are in
`outputs/hops_pulse_release_20261025_v1/`. The new physical run is
`outputs/lab_runs/hops_pulse_release_v1/`.

That run completed 5 µs. It retains 581/600 ions (96.83%), but ion RMS radius
grows 3.8576→6.1973 cm, median radius 3.8821→4.6147 cm, and core population
(r<2 cm) falls 39→26. It does not show repeated coherent compression. After
release, electrons fall from 7,847 at 0.2 µs to 7 at 0.3 µs. Electron retention
is a mechanism diagnostic, not the goal: the goals are ion density, inward
motion, useful ion energy, retention and measured input work.

The next matched experiment, `outputs/lab_runs/hops_pulse_source_on_v1/`, keeps
the same finite 0.2 µs pulse and ideal 1 mA electron source through 1 µs.
It isolates the source cutoff; ongoing source power remains in the ledger.
It is not free relaxation or a selected HOPS frequency. Frozen source/config
and launch are in `outputs/hops_source_control_20261027_v1/`. Compare ion
phase space/core density and potential after release before a timing scan.

The viewer now derives native E lines, species J/density/charge maps, ion radial
phase space and population history from all saved particles. Optional B is
only a softened plasma-current contribution, excluding feed/return/electrode/
displacement currents, steel magnetism and magnetic forces. Read-only evidence
and independent linear-E/current-loop comparisons are in
`outputs/hops_view_delivery_20261027/`. No physical solver change is implied.

Finer playback rerun: `outputs/lab_runs/hops_pulse_release_1ns_v1/` completed the
same 0–5 µs experiment with 1 ns full-state archives, ten times finer than the
original 10 ns cadence. Integration remains 5 ps; the old viewer skipped 2000
solver steps between frames, versus 200 now. No interpolated states or physical
parameter change. All native arrays match bitwise at all 501 common original
timestamps, including 5 µs; all 5001 dense states are retained. Source/config/launch evidence is in
`outputs/hops_temporal_archive_20261027_v1/`. This measures archival fidelity,
not temporal/spatial/particle convergence of the physical prediction.

**Measured pulse response:** `outputs/hops_relaxation_20261101_v1/` contains a
full native-ion measurement CSV, figure and JSON from `examples/hops_relaxation.py`.
After the 0.2 µs off command, RMS radius rises 10% in 1.315 µs. The population
inside the wire-centreline convex hull falls 10% in 2.401 µs and 20% in
3.923 µs. These first-sample times have 1 ns archive brackets, not converged
physical uncertainty. Hull faces are open and diagnostic only. At 5 µs,
446 ions remain inside that hull, 135 live outside it, and 19 are absorbed
(11 cage, 8 chamber). Of the 135 outside, 116 are beyond the entire 7.5 cm
outer cage radius. Thus 96.83% chamber survival is not 96.83% cage occupancy.

Instantaneous RMS dilation stays positive throughout the 4.8 µs post-release
record, and ends at 7.954 km/s. The final R/Rdot is 7.79 µs, an instantaneous
expansion scale, not a damping time. Neither half-empty nor e-fold cage
occupancy is reached. No bounded post-release breathing cycle, period or
single relaxation constant is identified. The fitted voltage fall is 20 ns,
reaching 1% of release voltage after 93 ns; this is the imposed boundary
response, not the ion response. Ion mean energy ends at 0.849 eV; a +3 kV
electrode is not a 3 keV ion-energy result. At zero cage voltage the centre
still has +0.843 V space-charge potential relative to the wall, which alone
does not demonstrate a restoring ion well. There is no collision operator
from which to infer a collisional/thermal relaxation time. Next extend the
source-free observation to actual cage half-empty times, tracking occupancy,
core density and phase alongside losses; do not choose an adaptive frequency
from this expanding transient alone. A concrete 0–25 µs single-pulse and matched
zero-bias/source-control design is now prepared in
[HOPS_RELAXATION_RUN.md](HOPS_RELAXATION_RUN.md), with executable configs and
waveform validation. Production cases are not launched. Force/particle/mesh
convergence remains open.

Two actual 20 ns +3 kV pulse/fall fixtures at 10/5 ps complete with finite
states. A 5 ns command with 2 ns response gives maximum particle-throughput
energy residual 7.257e-6→1.866e-6, charge error below 9.70e-27 C; emitted
charge/energy remain constant after source cutoff. Actual release is delayed
by 10/5 ps through fixed-mode phase arithmetic. Relative to the recorded
command time, voltage matches the analytic exponential within 2.8e-12 V.
The stricter comparison to an unquantized edge is retained separately (9.12→
4.55 V error). Poisson emission draws differ with dt, so these fixtures establish
budget/waveform behavior rather than strict particle phase convergence.

Measure radius, radial velocity, principal radii, potential, retained charge
and full energy/loss ledgers through release. A turning point can be biased by
particle loss; it must not be interpreted alone as coherent compression or
used to invent a relaxation frequency. Future re-pulse timing should use the
measured response rather than wait indefinitely for a prescribed contraction.

The original worker stopped at 2.52 µs on Oct 23, with unknown cause. On Oct 25,
its source RNG was recovered from the exact initial draw and emission stream,
matching birth counts, charge and injected energy at all 252 saved intervals.
An actual full-solver 2.51→2.52 µs replay then matched every saved particle and
field array bitwise, all physical diagnostics exactly, and the independently
reconstructed RNG/controller state. This establishes restart fidelity, not
physical convergence. Evidence: `outputs/hops_resume_20261025_v1/`, especially
`reconstruction.json` and `replay_verification.json`.

After recovery, v1 advanced to 2.628245 µs and was stopped at a full-step
checkpoint for the user's pulse/relaxation priority. All 264 states and original
source remain; the old DC/fixed tail will not start. Complete checkpoints retain
particles, fields, source RNG, controller and accumulated budgets. No feedback
switch occurred. The release radius was 3.0896 cm, but the filtered cloud radius
never fell below approximately its 3.8621 cm warmup reference through 2.59 µs.
This run therefore measured a voltage hold rather than a pulse response.
