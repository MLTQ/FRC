# crosscheck_bank.py

Independent lumped-circuit check of BUILD_DOSSIER section 2. It imports no
project circuit implementation. `python -m examples.crosscheck_bank OUT`
writes `OUT/main_bank.json` and prints the four comparison rows.

- `BankCase` holds the stated coil inductance, vacuum-current target and
  dossier values. The -0.05 T bias current is one quarter of the 0.2 T target;
  this checks the equivalent circuit, not the coil field calibration.
- `integrate` solves a bias branch, eight equal main modules represented in
  parallel, and a crowbar branch after the first main-capacitor voltage zero.
  L = diag(Lbranch) + Lcoil 11ᵀ, L I' = V−RI, V' = −I/C. At the first
  subsequent main-current zero the one-way start switch opens.
- `calculate` records peak current, slew, switching times, reversal and the
  energy identity including resistor dissipation. The ideal pre-crowbar
  quarter period is an independent limiting comparison.
- `main` finds the charge producing each stated vacuum-current target.

Assumptions: equal modules, main leads 30 nH/module, crowbars 8 nH and
2.1 mΩ/module; bias 2000 µF, 2.4 µH, 8 mΩ, initially at its current maximum.
Main lead/switch resistance is zero. Plasma loading and hardware design are
outside this calculation. Returned peaks are evaluated on dense saved ODE
solutions, not inferred from bank energy alone.

Measured four-case comparison: every dossier difference <2.1%; energy
identity residual <4.8e-11 J. The ideal quarter-period estimate agrees within
0.015%. Full values and limits are in `docs/CROSSCHECK.md`.
