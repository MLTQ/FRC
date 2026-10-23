# partial_pi_formation.py

Step 3 of bead frc-ibl.4. It asks how a partially pre-ionized fill changes FRC
formation in the benchmarked 2D formation path (coupled banks, exact vacuum
outside the tube), using the static neutral reservoir in
`plasma_bench/neutral_reservoir.md`. The module docstring gives the fill split
and the control.

## Runs (option A, `frc12_uniform_h2`: H2 3.4 mTorr, 0.2 T, 24×144, 0–12 µs)
- **Fill split:** the plasma starts as the fraction α of the fill under the
  coil; the rest of the fill under the coil is the neutral reservoir.
- **Ionization cost:** E_c = 30 eV per ion (about 2× the thermochemical
  minimum), with 60 eV as a sensitivity.
- **Control:** α = 1. It agrees with the Oct-20 design-scan run
  `d12_uniform_h2_p020` to 0.1% (flux 2.294 against 2.296 mWb at 1 µs). That
  run used an older solver revision, so the α = 1 run is the control here.
- **Bookkeeping, every run:** fluid energy residual ≤ 5e-11 J; response
  residual 0.12–0.13 J (1e-4 of the source work, the same as the control);
  mass change ≤ 2e-13.

| α (E_c) | neutrals left at 1 / 2 / 4 µs | ionization cost (J) | trapped flux at 4 / 12 µs (mWb) | separatrix at 4 / 12 µs (m) | length at 4 / 6 µs (m) | density at 4 µs (m⁻³) | Te+Ti at 2 / 12 µs (eV) | source work (J) |
|---|---|---|---|---|---|---|---|---|
| 1 (control) | – | 0 | 2.124 / 1.938 | 0.102 / 0.108 | 1.44 / 0.92 | 6.4e20 | 122 / 158 | 1309 |
| 0.3 (30 eV) | 26% / 6% / 3% | 64 | 2.143 / 1.862 | 0.094 / 0.105 | 0.83 / 0.55 | 8.2e20 | 97 / 155 | 1224 |
| 0.1 (30 eV) | 29% / 8% / 2% | 83 | 2.060 / 1.747 | 0.094 / 0.101 | 0.73 / 0.52 | 8.3e20 | 92 / 155 | 1164 |
| 0.03 (30 eV) | 31% / 10% / 2% | 89 | 1.941 / 1.605 | 0.091 / 0.096 | 0.76 / 0.52 | 8.3e20 | 90 / 148 | 1109 |
| 0.1 (60 eV) | 32% / 10% / 2% | 166 | 1.899 / 1.607 | 0.085 / 0.099 | 0.73 / 0.49 | 9.2e20 | 79 / 144 | 1127 |

## Reading
- **A closed FRC forms by 1 µs at every ionized fraction down to 3%.**
- **The neutrals do not survive the implosion.**
  - The hot sheath (60–120 eV) ionizes about 70% of the reservoir within
    1 µs, about 90% by 2 µs and all of it by 8 µs.
  - Nearly the whole fill ends up in the plasma: final mass is within 1% of the
    control.
- **The penalty is flux and size, not particle inventory.**
  - By 12 µs the trapped flux is 4% (α = 0.3), 10% (0.1) and 17% (0.03, or 0.1
    at 60 eV) below the control.
  - The FRC is smaller (separatrix −3 to −11%), denser, and much shorter during
    formation: 0.5–0.8 m against 0.9–1.4 m at 4–6 µs.
  - It is cooler early (Te+Ti −20 to −35% at 2 µs) but nearly recovers by
    12 µs.
- **The ionization energy is small: 64–166 J, 2–5% of the 3.2 kJ main bank.**
  - It is taken from the forming plasma when it is coldest.
  - Ions born at rest load the imploding sheet, converting 74–88 J of flow
    energy to heat.
  - The banks also deliver 6–15% less work.
- **These penalties are optimistic.**
  - No charge exchange: a real fill loses momentum and energy to fast neutrals
    that escape to the wall.
  - No neutral flow.
  - The rate is per atom-equivalent from H2 BEB.
  - E_c is declared.

## Design implication (frc-ibl.4)
The dossier's 60 J PI probably reaches α of order 0.1–0.3 (`pi_budget.md`).
- In this model that still forms a closed FRC, with roughly 4–10% less trapped
  flux and a shorter FRC than full pre-ionization.
- **A larger PI bank (about 0.7 kJ) is worth having but is not required for
  first formation.**
- Start with the PI as built. Upgrade if the bench shows a low density at main
  fire or low trapped flux.
- Charge exchange, which this model omits, would make the case for the larger
  PI stronger.
