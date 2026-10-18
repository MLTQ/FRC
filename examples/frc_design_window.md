# frc_design_window.py

Runs the first-cut FRTP design window for the user's tube:
`python -m examples.frc_design_window [OUTPUT]`. The default output is
`outputs/frc_design_window`. It takes a few seconds.

It writes:
- `checks.json`: independent checks of the rigid-rotor flux/decay closed forms
  (quadrature, finite difference) and of the sheet-coil field/inductance (sum of
  analytic elliptic-integral loops);
- `scan.csv`: H2/D2/N2 × fill 3–40 mTorr × B_e 0.1–0.75 T, flagged against
  `WINDOW`: 1 ≤ s_bar ≤ 4, T_e ≥ 10 eV, bank ≤ 5 kJ, ≤ 25 kV;
- `sensitivity.csv`: each empirical assumption at its range ends, stray
  inductance 15/60 nH, and 5/10 cm bores around the nominal point
  (D2, 10 mTorr, 0.3 T);
- `drive_options.csv`: series turns × parallel modules at the nominal fill and
  field for H2 and D2 (current per switch, charge voltage, capacitance per
  module, total energy, loop voltage per turn);
- `large_bore.csv`: a 12" ID bore with a 1.2 m coil and eight modules, H2/D2 at
  s_bar 2.5 and 3.5, for B_e 0.2–0.5 T;
- `summary.json`: device, assumptions, window and ranges.

The window limits are design preferences, not physics thresholds. Results and
their interpretation are in `docs/FRC_DESIGN_WINDOW.md`.
