# dual_energy_verification.py

`python -m examples.dual_energy_verification [OUTPUT]` (about 10 s; 16 cells).
For each configuration (off, default, no detector, detector 0.05, reset mode,
threshold 1e-3, bounded traces, dual energy plus bounded traces) it runs three
fixtures:
- **Standing cylindrical magnetic shock** (`standing_shock`, to t = 0.02). The
  independent reference is the exact steady jump, using `shock_errors` from
  `conservative_measurements`.
- **Coil pulse** (to 0.1). The final state is compared with the off run; a
  warm flow should be bit-identical.
- **Low-β coil pulse** (gas pressure 1e-6, β ≈ 2e-6, no physical or interface
  resistive heating, to 0.1). The independent reference is the adiabatic
  invariant A = p/ρ^γ, which stays at its uniform initial value in smooth flow.

Results go to `outputs/dual_energy_verification/results.json`.

Measured 2026-10-19, default configuration:
- shock pressure L1 2.5866 (off: 2.5853);
- coil pulse identical to off;
- low β completes, invariant error 3.0%, fluid energy residual 1e-16 J;
- off fails at t = 0.004.

Separately, the off path reproduced the pre-change commit bit for bit on both
fixtures.
