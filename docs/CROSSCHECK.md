# Independent build-dossier checks

## Main bank

The independent three-branch calculation in `examples/crosscheck_bank.py`
checks the equivalent circuit without importing `plasma_bench/coil_circuit`.
It uses the dossier's 77/93 nH coil inductances, eight 6 µF modules, 30 nH
leads per module, and the stated vacuum-current targets. The initial bias is
−0.05 T, at its current maximum; field per current is not independently
verified here (frc-ibl.15.3).

Each active branch obeys L I' = V−RI with
L = diag(Lbranch) + Lcoil 11ᵀ. Capacitor charge evolves with C V' = −I.
Crowbars close at the first main-capacitor voltage zero, and the unidirectional
main switch opens at the next main-current zero. Inductive/capacitive energy
plus integrated I²R is independently measured across both switch intervals.

Reproduce with `python -m examples.crosscheck_bank outputs/crosscheck_bank`.
Measured results are retained in `outputs/crosscheck_bank/main_bank.json`.

Values below are **independent / dossier**. The voltage is solved for the
stated vacuum coil-current target; this is not a plasma-loaded formation run.

| Case | Charge (kV) | Peak/switch (kA) | Initial dI/dt (kA/µs) | Crowbar (µs) | Reversal magnitude (kV) |
|---|---:|---:|---:|---:|---:|
| A, 0.2 T | 11.518 / 11.6 | 35.625 / 36 | 18.375 / 18 | 3.0467 / 3.05 | 2.847 / 2.9 |
| A, 0.3 T | 16.147 / 16.2 | 49.942 / 50 | 25.759 / 26 | 3.0466 / 3.05 | 3.991 / 4.0 |
| B, 0.2 T | 15.224 / 15.2 | 43.151 / 43 | 20.401 / 20 | 3.3243 / 3.33 | 3.449 / 3.5 |
| B, 0.3 T | 21.314 / 21.3 | 60.416 / 60 | 28.561 / 29 | 3.3242 / 3.33 | 4.829 / 4.8 |

All discrepancies are below 2.1%; none exceeds the 3% review threshold.
The largest difference is the rounded initial slew (18 versus 18.375).
Bank energies are 3.184, 6.257, 5.562 and 10.903 kJ. The main switch opens
about 0.755 µs after crowbar, leaving the negative capacitor voltage in the
table; zero capacitor voltage is not the switch-opening event.

With a slowly varying bias branch the independent quarter-period estimate is
π/2 sqrt[Cmain (Llead + Lcoil || Lbias)]: 3.04634 µs for A and 3.32381 µs
for B. The full bias/crowbar integrations differ by less than 0.015% before
crowbar. Their worst energy residual is 4.8e-11 J, measured across both
switch intervals. These checks validate the stated lumped equivalent numbers,
not field calibration, plasma loading, component ratings or a hardware design.
