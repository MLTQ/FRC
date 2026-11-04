# Goal status

Single current-state page. Keep it short: replace sections instead of appending
dated logs. The full pre-2026-10-18 log (about 6,000 lines), NEXT_SESSION and
THREAD_HANDOFF are preserved at git tag `archive/pre-refocus-2026-10-18`
(`git show archive/pre-refocus-2026-10-18:docs/GOAL_STATUS.md`).

## Goal (user, 2026-10-18)
Use the simulator to specify a real benchtop FRC experiment for the user's lab:
the fields (where and when), fill, pre-ionization and timing that form an FRC
torus, preferring longer lifetime. The user builds the pulsed power (budget
under 10k USD; hand-built capacitor banks and coils are fine) and has
pulsed-power experience and a radiation-safety consultant. The repository
specifies the required field and current waveforms and the equivalent
coil/bank numbers; it does not engineer the hardware.

Lab envelope:
- Glass tube about 15 cm diameter × 1 m. The bore is still open. McMaster's
  borosilicate tops out at 5.5" OD × 7" long (checked 2026-10-18), so the tube
  must come from a scientific-glass or process-pipe supplier. Bore options are
  in FRC_DESIGN_WINDOW.md.
- Two-stage roughing pump (~10 mTorr base); turbo planned.
- Gases: H2 as workhorse, D2 sparingly, N2 for cheap commissioning demos.
- Diagnostics: interferometer, borrowed high-speed camera, home-built magnetic
  probes and flux loops, and pressure, voltage and bias gauges.

## Plan (tracked under beads epic frc-ibl)
1. **Done: first-cut design window (frc-ibl.1)**, in
   [FRC_DESIGN_WINDOW.md](FRC_DESIGN_WINDOW.md). Nominal point: H2 or D2 at
   10 mTorr, external field 0.3 T, bias about −60 mT, main rise about 1.2–1.7 µs
   from a 16 cm × 50 cm single-turn coil (125 kA-turns). The coil is driven by
   four parallel modules, each about 2.7 µF at 9 kV carrying 31 kA (~400 J in
   total), with a crowbar. The estimates are T_e+T_i ≈ 90 eV,
   s̄ ≈ 2–2.8, r_s ≈ 3.5 cm and a classical τ_φ of at most ~18 µs. The
   empirical formation efficiencies are the main uncertainty.
2. **Done: benchmark (frc-ibl.2)**, in [NUCTE_BENCHMARK.md](NUCTE_BENCHMARK.md).
   - NUCTE-III case 6 is reproduced with a coil matching its mirror ratio and
     field index, driven by its banks through the circuit model with plasma
     loading: r_s 8.8 cm, l_s 41–44 cm, Te+Ti 100–105 eV and flux 1.1 mWb at
     25–30 µs, against the published 8.0 cm, 38 cm, 120 eV and 1.1 mWb.
     Density is 35–40% low.
   - Formation converges between 16×128 and 32×256.
   - The 2D model does not predict lifetime, which is set by the n=2 mode.
3. **Done: coil circuit (frc-ibl.5)**, in `plasma_bench/coil_circuit.py` and
   FRC_DESIGN_WINDOW.md.
   - **Model:** equal-flux coil inductances (NUCTE coil 68.7 nH against the
     published 70 nH), banks, unidirectional switches, crowbars, and
     self-consistent plasma loading.
   - **12" bank:** 0.3 T from a −0.05 T bias needs 16.2 kV and 6.3 kJ on
     8 × 6 µF. 0.2 T needs 11.6 kV and 3.2 kJ.
   - **Crowbar:** ≤ 8–20 nH and about 2 mΩ per module.
   - **Timing:** a common collector tolerates 100–200 ns of skew.
   - **Pre-ionization bank:** about 10 µF at 3.5 kV, ringing at 150 kHz.
   - **Quadrupole for n=2:** about 7% of B_e, i.e. 8–15 kA in four bars
     (literature-based).
4. **First pass done: design scan (frc-ibl.6)**, in [DESIGN_SCAN.md](DESIGN_SCAN.md).
   Eleven coupled-bank formation runs in the 12" geometry.
   - **Option A (start here):** uniform coil, H2 at 3.4 mTorr, 11.6 kV
     (3.2 kJ) for 0.2 T. A closed FRC forms by 2–6 µs: 11 cm × 0.8–1.1 m,
     150 eV, 2 mWb.
   - **Option B:** mild graded coil (20 cm centre, 16.5 cm ends), D2. At 0.2 T,
     15.2 kV and 5.6 kJ give 11.5 cm × 0.6–0.7 m, 220 eV, 2.2–2.5 mWb.
   - **D2 on a uniform coil** forms late or not at all.
   - **Robustness:** results are insensitive to the gap floor and to
     pre-ionization temperature (1–5 eV).
   - **Expected diagnostic signals** are tabulated in the scan doc.
   - **Coil shape:** the user chose to keep both options as presets
     (2026-10-19). A 36×216 check of option A agrees within 10–15%.
5. **Experiment playground (frc-ibl.9)**, in
   [EXPERIMENT_PLAYGROUND.md](EXPERIMENT_PLAYGROUND.md).
   - **Presets:** named experiments (`preset=NAME`), extra banks with their own
     firing times, aux coils, and a coaxial insulating inner tube with coils
     inside it.
   - **Inner-tube solver:** verified against an independent 1D diffusion
     solution at second order.
   - **3D viewer and comparison page:** `examples/formation_display.py`.
   - **Demo results:**
     - A re-boost compresses but adds no flux.
     - End plugs barely matter, because the crowbarred main coil conserves its
       flux.
     - A 3" inner tube costs ~40% of the flux, and the separatrix ends up on
       the tube.
     - A slow flux core inside the tube grows φ to 2.8 × the control by 40 µs,
       but makes the FRC short and fat.
     - A fast core punch splits the FRC.
   - **Lifetime levers (2D, to 40 µs):**
     - 0.3 T collapses after 20 µs, from extra open-line conduction loss.
     - A slow bias-time mirror forms no real mirror and adds no flux.
     - A main-coil re-boost restores the field but not the FRC's shape.
     - A 1 kV flux core (300 µF, 150 J, fired at 4 µs) holds φ flat at
       1.33–1.41 mWb, but the FRC still shortens to E ≈ 1.2, because it keeps
       losing energy.
     - Plain option A at 0.2 T remains the lifetime baseline: 1.11 mWb and
       E 2.6 at 40 µs.
     - Tilt and n = 2 are not modelled.
6. **Build dossier:** [BUILD_DOSSIER.md](BUILD_DOSSIER.md) gathers every coil's
   voltage, current, energy and timing, plus the pre-ionization, bias,
   quadrupole and shot sequence.
7. **3D test-particle orbits (frc-ibl.12)**, in EXPERIMENT_PLAYGROUND.md
   ("Particle orbits") and `examples/particle_tracking.md`.
   - **Method:** relativistic Boris full orbits in the saved fields plus a 3D
     quadrupole and RMF antenna, verified analytically; p_θ converges at
     second order.
   - **Kinetic ions:** they are kinetic (≈60% cross the null within 2 µs).
   - **Option A:** keeps 95% of its ions from 10 to 30 µs.
   - **Quadrupole at 7% B_e:** loses an extra 24–49% to the wall in m=2
     sectors at |z| ≈ 0.5–1.1 m. This is a vacuum-superposition upper bound.
   - **Multipole scan (option A):** the quadrupole costs 7/15/24/31 points
     at 3/5/7/9% of B_e. An 8-bar octupole at 7% costs 4.8 points; its n=2
     threshold needs frc-ibl.11.4.
   - **3" inner tube:** absorbs 13–23% of ions promptly.
   - **RMF (1 MHz, 5 mT):** entrains electrons near the null at 0.93 ω, but
     at 0.2 T an RMF needs f ≳ 4 MHz to exceed ω_ci.
8. **Remaining gaps:**
   - pre-ionization ionization fraction (frc-ibl.4); temperature is now known
     not to matter;
   - N2 commissioning expectations (frc-ibl.7);
   - reset corrections at low fill (frc-ibl.2.10);
   - Hall-cutoff sensitivity and a resolution check of the Hall run
     (frc-ibl.13).

**Target geometry (user, 2026-10-18): 12" bore to maximize lifetime.**
- Tube about 2.2 m long; 1.2 m single-turn coil fed by eight modules.
- Baseline fill about 3.5 mTorr. The design scan now favours H2 on a uniform
  coil, or D2 with a graded coil.
- Bank charge, superseding the earlier lossless estimate: first light at 0.2 T
  is 11.6 kV and 3.2 kJ; 0.3 T needs 16.2 kV and 6.3 kJ, both from a −0.05 T
  bias on 8 × 6 µF. In 2D, 0.3 T is hotter but shorter-lived.
- Classical τ_φ is 113–354 µs. LSX/FRX-C scaling suggests about 50–150 µs is
  realistic with quadrupole stabilization.
- See FRC_DESIGN_WINDOW.md. The 6" nominal point above is kept for comparison.

## Solver state (honest summary)
**The formation path** is `ConservativeVacuumMHD`, driven by
`examples/nucte_formation.py` with `drive="coupled"`:
- 2D axisymmetric single-fluid MHD with exact vacuum outside the tube;
- Spitzer + Chodura resistivity;
- energy-consistent dual energy and the compatible energy path;
- coil banks coupled self-consistently;
- optional 2.5D (v_θ, B_θ), a declared Hall term with a density cutoff, and a
  threading-current boundary for a rod in the inner tube (frc-ibl.10,
  2026-10-20). Off by default and bit-identical when off.

It is benchmarked against NUCTE-III for formation, and formation converges
between 16×128 and 32×256.

**Hall result (option A, frc-ibl.10).** With Hall on, the FRC forms at the
same time but contracts axially 0.6–1.8 µs sooner, to the same 64 cm length.
r_s stays within 3% and trapped flux within 4% of the 2D run. B_θ reaches
85 mT, and J_z reaches 65% of peak J_θ. The rebound near 15 µs is larger (see
NUCTE_BENCHMARK gap 5). The cost is about 2.2× the 2D wall time.

**3D hybrid work (frc-ibl.11).** All six requested Cartesian hybrid shots
completed actual 0–40 µs from initial fill: baseline, end plugs, re-boost,
inner tube, flux core and core punch. Each has 40,000 ion steps and 401
audited viewer frames. The particle slider selects 0–4,000 saved dots only.
There is no azimuthal filter, mirrored sector, artificial mode or energy/
temperature rescaling. Initial fill and coil geometry are specified inputs;
the gold n=0 outline is an averaged diagnostic, not the evolved 3D topology.

**Higher-resolution reruns (frc-ibl.11.10).** The first 10 mm CPU baseline
failed at 1.141 µs because the field timestep did not bound outgoing electron
energy transport at cut walls. Its queue was stopped with evidence retained.
The shared CPU/Metal energy operator now takes positive conservative internal
substeps; independent analytic translation conserves energy to roundoff and
converges. No energy/temperature rescale is used. The corrected six-shot queue,
`outputs/hybrid3d_shots_r10_cfl_metal/`, and one matched CPU64 baseline,
`outputs/hybrid3d_shots_r10_cfl_cpu_control/`, use immutable
sources: 10×10×10 mm, eight ions/cell, 0.5 ns and actual 0–40 µs. Baseline has
1,290,080 ions (16.3× coarse). The GPU queue holds other cases until baseline
passes 1.3 µs; then at most two shots use two CPU threads each. The CPU control
had one worker. Both corrected GPU/CPU baselines completed 40 µs. An Oct 25
audit finds the GPU coordinator and both remaining workers absent: end plugs
stopped at 17.2 µs and re-boost at 5.2 µs on Oct 22, with no traceback explaining
the exit. Their saved endpoints are finite. Inner tube, flux core and core punch
never started and have no active queue. Stale running labels are corrected;
evidence is in `outputs/hybrid3d_status_20261025/interruption.json`. These FRC
snapshots save 4,000 sampled ions, not complete restart state; exact continuation
from those frames is unavailable. Preserve all outputs before any new replay.
`DONE` means terminal queue states. The new energy timestep means this is not
a pure grid comparison with v4. Numerical-loss and precision limits remain open.

**GPU prototype (frc-ibl.11.11).** Optional Metal transport and CIC deposition
are built. At 1.29 million ions, deposition is 3.18× faster including transfers;
complete startup steps take 14.42 s versus CPU 16.57 s (1.149× throughput).
Independent comparisons converge, and a coupled 0.75 µs replay preserves the
load budget within 6.1 µJ of CPU. The smaller-grid replay has essentially
equal runtime. Core integration/vacuum remain CPU64; GPU blocks explicitly
use float32. The fine CPU/GPU control shows response-face L2 differences
of 23% at 3.9 µs despite close global budgets; above-cutoff nodal field and
density differences have both reached 16%. Particle-count and physical-loss
convergence remain unresolved. Partial saved-frame previews of the fine
baseline, CPU control and end plugs are in the viewer's Run menu with actual
time spans; they do not claim completion or refresh automatically.
Full-shot GPU accuracy/performance are not established. Claude accepted the
physics lane on frc-ibl.11.9; Codex retains performance/run coordination.
Both fine baselines completed 40 µs, with numerical removal 982/969 J against
1,468/1,469 J source work (GPU/CPU). Consult the queue manifest for remaining
configurations. At matched 10.8 µs, trapped flux is 1.377 versus
0.901 mWb and length 0.72 versus0.54 m for GPU/CPU despite nearly identical
numerical removal (~746 J). At 11.6 µs baseline removal reaches 762 J, 51% of
source work; energy-ledger error is 0.44 J. Numerical loss and divergence between
trajectories prevent configuration rankings. Evidence: `outputs/hybrid3d_status_20261022/`.
See HYBRID_3D.md for evidence and limits.

**Backend isolation (frc-ibl.11.10.1).** Five actual fine-grid 0–0.1 µs controls
isolate the early difference: CPU repeat is bitwise identical, transport-only
response difference 7.16e-11, deposition-only/combined 4.869e-6. A separate
CPU64-cell/Metal32-local-fraction CIC prototype completes two matched coupled
runs, reducing response difference to 5.61e-9 (~870x); pressure 9.84e-9 and density
1.08e-7. Independent number/linear-field-work comparisons were measured; the
work error and numerical cooling are not cured. Normal production imports and
all running campaigns keep their original source. The frozen combined prototype
now completes 1 µs: response difference is 1.419% versus the original GPU's
9.501% (6.70x reduction), with dense-region nodal field difference 0.0109%.
Low-density differences still grow; numerical removal remains about 157 J
against 1,002 J source work. CPU/original reference archives are float32.
Next: a same-seed corrected replay to measure residual atomic sensitivity
before longer adoption. Artifacts: `outputs/hybrid3d_precision_*` and
`outputs/hybrid3d_cic_coordinates_endpoint_v2/`. Late-time/physical convergence
remains open; Claude owns the physical-loss/interface lane.

**CPU array optimization (frc-ibl.11.12).** Opt-in float64 fused stencils
use about 19% less time in the measured fine GPU replay; CPU fine uses
19–31% less under varying host load. Complete CPU64 fields, particles, heat
and flux arrays match the original exactly at both measured endpoints, and
independent comparisons retain their errors/conservation. New launches can
select `arrays=1` or campaign `--arrays`; current queues keep their saved
source and progress. This speeds evaluation without curing physical losses.

**Low-density diagnosis (Claude, frc-ibl.11.9).** Snapshot post-processing
locates numerical Joule dissipation in under-resolved/carrier-free gas nodes;
the deletion keeps them cold and resistive. His separate 0–1 µs vacuum-mask
prototype reduced the Joule-channel sink but deleted more field energy during
mask projection: total numerical removal rose from 68 to 112 J. It is not
adopted. His fixed-mesh 1 µs comparison at actual counts 8→27→64 reduces
numerical heat removal 68→53→45 J, with close n=0 geometry. Full 3D convergence
does not follow: above-cutoff field differences remain 12% for 27→64 and
thin-node differences 51%. The old loader silently rounded nominal 16/32 to
27; new launches require exact cubic counts and record their lattice. His
passive carrier ledger records orphaned-edge Ohmic heat 43 of 95 J in
startup and 111 of 113 J in analytic compression. Blocking electron-energy
transport into empty nodes alone gives startup removal 74 versus 68 J: not
adopted. His completed 5/2/1% floor comparison reduces removal 68→48→39 J,
but raises empty-node removal and boundary enthalpy outflow. The early
midplane radius changes 0.121→0.112→0.112 m; this single scalar at 1 µs does
not establish full 3D convergence or justify changing the running campaigns.
Next bounded comparison pairs the 1% floor with the electron-edge block,
retaining the original controls and all energy channels. This is a numerical
closure study, not a physical-vacuum solution. See
`plasma_bench/hybrid3d_carrier_ledger.md`,
`examples/lowdensity_loss_check.md` and `examples/vacuum_mask_replay.md`.

The full-shot path includes exact exterior/inner vacuum, absorbing walls,
banks coupled to full 3D current and conservative electron internal energy
with Spitzer heating. Independent contact, heating, loop reciprocity and
boundary-mode comparisons pass or converge; prior energy-failing campaigns
are marked INVALIDATED. See [HYBRID_3D.md](HYBRID_3D.md) for measured results.
The v4 accounting residual is at most 0.46 J after the switched-bank flux
reference is included, but explicit numerical heat removal totals 527–780 J
(43–60% of load work). Early combined grid/time refinement changes magnetic
and ion energies substantially. These are development runs, not converged
formation, stability, confinement or lifetime predictions. Next: measure the
corrected high-resolution comparison and resolve density-cutoff/empty-node losses before
interpreting configuration rankings; M3/M4 and the overarching goal stay open.

**The 2D formation path still lacks:**
- two-fluid physics beyond the Hall term: the ∇p_e (battery) term, electron
  inertia and two temperatures;
- neutrals and ionization, so pre-ionization is assumed;
- 3D modes (n=2), so it cannot predict lifetime.

**Other repository paths** (PROJECT_OVERVIEW.md):
- The finite-cylinder and Lagrangian paths are not used for formation.
- The argon slab and Vlasov work are paused (frc-transport.*, P3).

**Positive-grid electrostatic side experiment (frc-cac).** The user's reference
hardware is a 40 cm chamber and 15 cm open truncated-icosahedron cage of ~2 mm
wire. Historical air at 20 mTorr, 0 to +3 kV and 120 Hz are reference conditions, not
constraints or a selected resonance. A separate full 3D actual-mass kinetic
electron/ion Vlasov–Poisson path is now built: conforming gas mesh around all
wires, self-consistent space charge, absorbing steel, optional grounded-wall
electron source, actual electrode work and uncorrected energy/charge budgets.
No prescribed ion well, neutralizing atmosphere, symmetry filter or energy repair.
The wall is explicitly assumed spherical; the source is an ideal emitter, not
an air discharge. The FRC quasi-neutral hybrid is not used for this model.

The final source-driven 0–50 ns app run completed 101 frames with full 3D state
archives (`outputs/lab_runs/iec_20261021_234301_0adb/`). Independent Gauss-law
Poisson comparisons converge, actual cage refinement approaches the independent
boundary solution, and paired 20 ns DC controls improve with 20→10 ps steps.
However small-orbit force errors remain 12.6% endpoint/5.05% energy at 10 mm;
particle/mesh convergence, local wire loss accuracy, collisions and ion-period
runs are unresolved. No POPS, fusion or optimal-parameter result exists. Next:
resolve force accuracy (frc-cac.2.1), then self-consistent well/free-ion/driven
comparisons over actual ion periods. See [IEC_3D.md](IEC_3D.md).

**HOPS (frc-cac.2.2).** User-defined **Harmonically Oscillating Plasma Sphere**
uses causal measured-radius feedback on the positive cage voltage. Actual
+3 kV source-off startups at 10→5 ps reduce particle-throughput energy error
0.378→0.113%, with charge error below 1.87e-27 C. Exact saved-interval replay
recovered the interrupted first run (evidence: `outputs/hops_resume_20261025_v1/`).
The user now prioritizes pulsing/relaxation over voltage holds. The old controller
waits for 20% inward contraction before release; this never occurred, so its
boundary stayed +3 kV. It is checkpointed/stopped at 2.628245 µs with all
264 states retained; its DC/fixed tail is stopped as well. This was time evolution
under a held electrode, not a demonstrated pulse/relaxation experiment.

The single-pulse experiment completed actual 0–5 µs at 5 ps, with 501 saved
frames/full states (checked Oct 27). It commands +3 kV for 0.2 µs, with a 20 ns
voltage response, then source-free relaxation. The ideal 1 mA emitter
also stops at release; source charge/energy stay explicit. It reuses the same
seed, mesh and full 3D particle/Poisson equations. Pulse duration is a declared
probe, not a found optimum. Evidence and frozen source:
`outputs/hops_pulse_release_20261025_v1/`. Measure actual radius/velocity turning
points, well evolution and loss before adapting subsequent pulse times.
The completed run is `outputs/lab_runs/hops_pulse_release_v1/`.
It retains 581/600 ions, but RMS radius grows 3.86→6.20 cm and core-ion count
falls 39→26. A matched 0–1 µs finite-pulse control with continuing ideal electron
source completed as `hops_pulse_source_on_v1` (paired analysis pending), to isolate source
cutoff from timing. Ion compression/retention and input work are the targets;
electron loss is an explanation for field evolution, not the objective.
The finer `hops_pulse_release_1ns_v1` completed 0–5 µs with 5001 native states
at 1 ns cadence and unchanged 5 ps integration/physics. All native arrays match
bitwise at all 501 original common timestamps. This verifies archival independence,
not physical convergence. Native response evidence is in
`outputs/hops_relaxation_20261101_v1/`: after release, RMS radius grows 10% in
1.315 µs and cage-envelope occupancy falls 10% in 2.401 µs. At 5 µs, 446 ions
are inside the open cage envelope, 135 live outside and 19 absorbed; chamber
survival is not cage confinement. Dilation stays outward after release; no
bounded breathing cycle or damping time is identified. Half-empty/e-fold cage
occupancy is not reached within the 4.8 µs release record. A concrete 0–25 µs
single-pulse/zero-bias matched-source pair is designed, with unchanged 5 ps
integration and 1 ns archives; no second pulse occurs within the interval.
Configs and waveform/reference evidence are prepared; production is not launched.
See [HOPS_RELAXATION_RUN.md](HOPS_RELAXATION_RUN.md). Next execute that longer
source-free response and inspect density/phase before choosing adaptive timing.
Collisional relaxation is absent from this development model.
Two short pulse/fall fixtures at 10/5 ps complete with finite states and
small energy/charge residuals; actual edge timing resolves to one step.
No radius-frequency law, harmonic well or sphere is forced. Inputs, ideal source and unresolved
force/collision limitations are in [HOPS.md](HOPS.md). No resonance or fusion
result exists yet; the live app plots saved radius/voltage in microseconds.

**Unified local app (frc-8qy, first version complete).** Fieldwork at
`http://localhost:8792` shares experiment/run navigation and displays both the
new IEC setup/run/replay and the existing FRC scientific viewers. The IEC
particle slider is display-only; real inputs, model limits, potential and raw
energy budgets remain visible. App launch, completion, cancellation, playback
and both 2D/3D FRC adapters were exercised. FRC campaigns retain their existing
run tools and ongoing progress. See [EXPERIMENT_LAB.md](EXPERIMENT_LAB.md).
Oct 27: the stopped FRC static viewer on 8791 was restored as a detached
loopback server (PID 20440), and its saved 3D rendering was browser-verified.
Fieldwork on 8792 was healthy and retained. Evidence: `outputs/web_servers_20261027/`.
The Oct 27 native-state viewer adds E lines, full-population species J/density,
ion radial phase space and live counts. Optional B is explicitly a softened
plasma-only reconstruction, excluding circuit/steel/EM-force physics. App 8792
was replaced by PID 28294; no simulation process was signalled. Evidence:
`outputs/hops_view_delivery_20261027/`.
Playback now retains the last complete frame while raw/derived data load and
commits all panels in one animation frame. Backing canvases resize only when
needed; revised asset URLs avoid stale scripts. No solver or saved data change.
The overarching FRC and electrostatic research goals remain open.

## Repository state
- Simulation outputs are no longer tracked; they regenerate from `examples/`.
  On 2026-10-18 the user had the old local outputs (~17 GB) deleted. Only
  `outputs/frc_design_window` and the user's `outputs/frc_showcase/poster.png`
  remain. Earlier figures linked from older docs must be regenerated with their
  example drivers.
- Outputs were purged from all history on 2026-10-18 at the user's request.
  The purge used git-filter-repo, which rewrites every commit hash; `.git` went
  from 10 GB to 9 MB. Codex's two turn-diff checkpoint refs were deleted, and
  the index of its worktree (/private/tmp/frc-finite-drift, branch
  codex/finite-drift-collisions, unmerged) was reset to the rewritten branch.
- Current simulation evidence is stored locally in `outputs/`. The corrected
  20 kA rod-before run completed 10 µs. All six hybrid v4 runs completed
  40 µs; `outputs/hybrid3d_full_shots_v4/` retains exact source archives,
  campaign/energy summaries and complete viewer audits. All eight older 2D
  scans resumed automatically after DONE; both priority-pause records
  (`/private/tmp/frc-3d-priority-pause-v4.json` and `-v4b.json`) say
  `resumed: true`. Invalidated v1–v3 data remain for diagnosis. Outputs are
  not committed, no remote is configured, and issue state is exported locally.
- Unrelated: a loopback http.server from an earlier agent may still serve
  `/private/tmp/frc-wrap-viewer-20261017` on port 8790.

## Rules in force
See [AGENTS.md](../AGENTS.md): simulation-first, no new unit tests, exact
vacuum outside the plasma, independent comparisons, and beads for tracking.
