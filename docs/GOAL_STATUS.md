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
2. **In progress: benchmark (frc-ibl.2).** See
   [NUCTE_BENCHMARK.md](NUCTE_BENCHMARK.md).
   - **Experiment and solver.** NUCTE-III (coil r = 0.17 m, D2 at 10 mTorr)
     on `ConservativeVacuumMHD`. The driver defaults are now the benchmark
     configuration: Spitzer + Chodura η, entropy-criterion reset dual energy,
     the compatible energy path, CT-Poynting magnetic energy flux and a 1%
     gap floor.
   - **Fixed in this pass:**
     - a dual-energy fault that created 4–5× the coil work as spurious heat
       (frc-ibl.2.8);
     - a response-energy error of up to 1.4× the coil work, removed by the
       compatible path (now ~1e-4; frc-ibl.2.9).
   - **What the model reproduces.** Case 6 (mirror 1.8) with a graded coil
     that matches NUCTE's R_M and field index, and drive values fixed from
     NUCTE's own numbers (bias −0.05 T, 0.3 T peak):
     - the sequence, with contraction to l_s 44 cm by 8 µs;
     - r_s 8.5 cm (published 8.0), l_s 44 cm (38), n̄ 2.6e21 (2.2e21) and
       T 162 eV (120) at 10 µs;
     - T 126–142 eV and n̄ 1.7–2.3e21 at 20–30 µs;
     - trapped flux 1.0 mWb (the paper's estimate is 1.1).
   - **Resolution:** formation quantities agree to about 10% between 16×128
     and 32×256.
   - **Lifetime is not predicted by this 2D model.** The 16×128 collapse at
     about 41 µs is numerical. At 32×256 the FRC's flux e-folds in about
     250 µs, while NUCTE's 50 µs lifetime is set by the n=2 rotational mode.
   - **Robustness:** the FRC is insensitive (≤13%) to the gap closure. The
     remaining counted energy corrections live in the gap (frc-ibl.2.10).
   - **Next action:** frc-ibl.2 is complete enough for formation design.
     Remaining gaps:
     - the circuit model (frc-ibl.5);
     - pre-ionization (frc-ibl.4);
     - an empirical/stability lifetime model, with quadrupole
       stabilization (frc-ibl.8).
     After those comes the design scan (frc-ibl.6).
3. Fill only the gaps the benchmark exposes:
   - hydrogenic plasma in the driven solver (frc-ibl.3);
   - pre-ionization and bias-flux trapping (frc-ibl.4);
   - the coil circuit with crowbar and end mirrors (frc-ibl.5).
4. Design scan and independent cross-check for the user's geometry (frc-ibl.6).
   WarpX is already installed on the user's GPU host; see ENGINE_COMPARISON.md.
5. N2 commissioning expectations (frc-ibl.7).
6. Quadrupole stabilization of the n=2 rotational mode (frc-ibl.8).

**Target geometry (user, 2026-10-18): 12" bore to maximize lifetime.**
- Tube about 2.2 m long; 1.2 m single-turn coil fed by eight modules.
- Baseline D2 at about 3.5 mTorr. First light at 0.2 T (8.3 kV, 1.6 kJ),
  upgrade to 0.3 T (12.5 kV, 3.7 kJ) on the same capacitors.
- Classical τ_φ is 113–354 µs. LSX/FRX-C scaling suggests about 50–150 µs is
  realistic with quadrupole stabilization.
- See FRC_DESIGN_WINDOW.md. The 6" nominal point above is kept for comparison.

## Solver state (honest summary)
- Several separate research paths exist, not one validated formation solver.
  See [PROJECT_OVERVIEW.md](PROJECT_OVERVIEW.md).
- Mapping which path can run a full FRTP formation is the first part of
  frc-ibl.2.
- The configured metre-scale coil experiments
  ([MATERIAL_EXPERIMENT.md](MATERIAL_EXPERIMENT.md)) used 0.1–20 mT fields on a
  4×16 mesh with argon. That is far below an FRC-forming regime, and the current
  was unconverged between grids (30–60% differences).
- The argon reaction/transport slab and Vlasov work are paused
  (beads frc-transport.*, P3). Resume only if the benchmark needs them.

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
