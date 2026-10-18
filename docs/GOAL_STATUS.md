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
- Pyrex tube (McMaster), about 15 cm diameter × 1 m.
- Two-stage roughing pump (~10 mTorr base); turbo planned.
- Gases: H2 as workhorse, D2 sparingly, N2 for cheap commissioning demos.
- Diagnostics: interferometer, borrowed high-speed camera, home-built magnetic
  probes and flux loops, and pressure, voltage and bias gauges.

## Plan (tracked under beads epic frc-ibl)
1. **Done: first-cut design window (frc-ibl.1)**, in
   [FRC_DESIGN_WINDOW.md](FRC_DESIGN_WINDOW.md). Nominal point: H2 or D2 at
   10 mTorr, external field 0.3 T, bias about −60 mT, main rise about 1.2–1.7 µs
   from a single-turn 16 cm × 50 cm coil (125 kA, ~600 J; ~7.5 µF at 12.5 kV
   with ≤30 nH stray), and a crowbar. The estimates are T_e+T_i ≈ 90 eV,
   s̄ ≈ 2–2.8, r_s ≈ 3.5 cm and a classical τ_φ of at most ~18 µs. The
   empirical formation efficiencies are the main uncertainty.
2. **Next: benchmark (frc-ibl.2).** Reproduce a documented small FRTP
   experiment with the most complete coil-driven solver path before trusting
   any design prediction.
3. Fill only the gaps the benchmark exposes:
   - hydrogenic plasma in the driven solver (frc-ibl.3);
   - pre-ionization and bias-flux trapping (frc-ibl.4);
   - the coil circuit with crowbar and end mirrors (frc-ibl.5).
4. Design scan and independent cross-check for the user's geometry (frc-ibl.6).
   WarpX is already installed on the user's GPU host; see ENGINE_COMPARISON.md.
5. N2 commissioning expectations (frc-ibl.7).

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
  The local `outputs/` directory still holds the old data (~17 GB).
- History before the untracking is unchanged, so `.git` is still ~10 GB.
  Purging outputs from history needs the user's go-ahead.
- No processes are left running. The paused 2058-step replay (PID 33081) was
  killed on 2026-10-18 at the user's request.
- Unrelated: a loopback http.server from an earlier agent may still serve
  `/private/tmp/frc-wrap-viewer-20261017` on port 8790.

## Rules in force
See [AGENTS.md](../AGENTS.md): simulation-first, no new unit tests, exact
vacuum outside the plasma, independent comparisons, and beads for tracking.
