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
6. **Build dossier:** [BUILD_DOSSIER.md](BUILD_DOSSIER.md) gathers every coil's
   voltage, current, energy and timing, plus the pre-ionization, bias,
   quadrupole and shot sequence.
7. **Remaining gaps:**
   - pre-ionization ionization fraction (frc-ibl.4); temperature is now known
     not to matter;
   - N2 commissioning expectations (frc-ibl.7);
   - reset corrections at low fill (frc-ibl.2.10).

**Target geometry (user, 2026-10-18): 12" bore to maximize lifetime.**
- Tube about 2.2 m long; 1.2 m single-turn coil fed by eight modules.
- Baseline fill about 3.5 mTorr. The design scan now favours H2 on a uniform
  coil, or D2 with a graded coil.
- Bank charge, superseding the earlier lossless estimate: first light at 0.2 T
  is 11.6 kV and 3.2 kJ; the upgrade to 0.3 T is 16.2 kV and 6.3 kJ, both from
  a −0.05 T bias on 8 × 6 µF.
- Classical τ_φ is 113–354 µs. LSX/FRX-C scaling suggests about 50–150 µs is
  realistic with quadrupole stabilization.
- See FRC_DESIGN_WINDOW.md. The 6" nominal point above is kept for comparison.

## Solver state (honest summary)
**The formation path** is `ConservativeVacuumMHD`, driven by
`examples/nucte_formation.py` with `drive="coupled"`:
- 2D axisymmetric single-fluid MHD with exact vacuum outside the tube;
- Spitzer + Chodura resistivity;
- energy-consistent dual energy and the compatible energy path;
- coil banks coupled self-consistently.

It is benchmarked against NUCTE-III for formation, and formation converges
between 16×128 and 32×256.

**What it lacks:**
- Hall/two-fluid physics, so the end reconnection rate in weak mirrors is
  uncertain;
- neutrals and ionization, so pre-ionization is assumed;
- 3D modes (n=2), so it cannot predict lifetime.

**Other repository paths** (PROJECT_OVERVIEW.md):
- The finite-cylinder and Lagrangian paths are not used for formation.
- The argon slab and Vlasov work are paused (frc-transport.*, P3).

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
- No processes are left running. The paused 2058-step replay (PID 33081) was
  killed on 2026-10-18 at the user's request.
- Unrelated: a loopback http.server from an earlier agent may still serve
  `/private/tmp/frc-wrap-viewer-20261017` on port 8790.

## Rules in force
See [AGENTS.md](../AGENTS.md): simulation-first, no new unit tests, exact
vacuum outside the plasma, independent comparisons, and beads for tracking.
