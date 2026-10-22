# hybrid3d_carrier_ledger.py

Bead frc-ibl.11.9, physical-loss/interface lane. This is the bounded comparison
agreed with Codex on 2026-10-22. It is a subclass in a new file; there are no
shared solver changes.

`CarrierLedgerHybrid` keeps `DeviceHybrid`'s field, Ohm's law and density floor
unchanged. There is no mask and no projection, so no current is deleted.

## Ledger (passive)
- **Ohmic heat by edge carrier class.** Every RK stage, with the solver's own
  (1, 1, 4)/6 weights, the Ohmic heat is booked by the macro-ions in each
  edge's four cells after the push:
  - orphaned: none;
  - sparse: fewer than `sparse_ions` = 4, under one per cell;
  - resolved: the rest.
- **Cutoff-removed heat by node class:** plasma (n ≥ n_floor), thin, or empty.
- **Edge current magnitude** Σ|J|·A by class, at diagnostic times.
- **Departures:** a count of cells whose last macro-ion left.
- **Checks:** the class sums are checked against the solver counters at every
  diagnostic row (to about 1e-13 J). With `electron_edge=0` the run is
  bit-identical to the baseline.

## Option `electron_edge` (default on)
- The energy closure becomes `EdgeEnergyElectrons`: no internal-energy flux
  crosses a face touching a node with no ions at the start and end of the ion
  step. Energy on such nodes is removed and counted.
- With no interface scale, its pressure force, pressure work and compression
  velocity are exactly the base max(n, n_floor). The only change is the
  transport block.

Results are in `examples/vacuum_mask_replay.md` (startup) and
`examples/edge_compression_check.md` (analytic compression).

## Result (2026-10-22)
- **The startup and analytic sinks are both Ohmic dissipation of carrier-free
  (orphaned) current:**
  - startup: 43 of 95 J of Ohmic heat by 1 µs, on edges that carry 35% of the
    current magnitude;
  - compression: 111 of 113 J.
- **Its dissipation is set by an undefined resistivity.** The density is zero
  there, and the closure's Te = p/(max(n, 1 m⁻³) e) gives either an accidental
  η ≈ 0 (electron energy leaked in during the step) or the 1 eV cap.
- **The electron-edge block removes the leak but exposes the cap.**
  - Startup: total removal 68 → 74 J.
  - Compression: 182 → 144 J.
- **It does not justify itself without a field-side change.**
- **Floor convergence** (5% / 2% / 1%, 0–1 µs, results in
  `examples/vacuum_mask_replay.md`):
  - Orphaned Ohmic heat falls 43 → 21 → 12 J, and total removal 68 → 48 → 40 J.
  - The 1 µs separatrix converges at 2% (0.112 m, against 0.121 m at 5%).
  - Empty-node removal (20 → 28 J) and wall enthalpy outflow (39 → 74 J) rise
    as the floor drops.

