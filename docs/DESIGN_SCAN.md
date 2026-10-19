# 12-inch design scan: formation runs (frc-ibl.6, started 2026-10-19)

## Purpose
Use the benchmarked formation model (NUCTE_BENCHMARK.md) to choose between
build options for the user's 12" device: gas, field, bias and coil shape. Then
state the build specification with its knobs and expected diagnostic signatures.

**What these runs predict:** formation, meaning reversal, end reconnection, axial
contraction, and the size, density, temperature and trapped flux of the
resulting FRC.

**What they do not predict:** lifetime. That is set by the n=2 mode and anomalous
transport (frc-ibl.8, FRC_DESIGN_WINDOW.md).

## Configuration
Driver: `examples/nucte_formation.py` with the benchmark defaults:
- Spitzer + Chodura resistivity;
- entropy-criterion reset dual energy;
- the compatible energy path and CT-Poynting magnetic energy flux;
- 1% gap floor;
- open-line conduction f = 0.1.

**Geometry:**
- Tube: 12" ID (r_w = 15.24 cm), 2.2 m long, reflecting ends.
- Coil: single turn, 1.2 m long, 24 straps of 5 cm on one collector:
  - uniform at r = 16.5 cm, or
  - graded from r0 at |z| ≤ 0.1 m down to 16.5 cm at the ends.

**Fill:** D2 at 3.5 mTorr or H2 at 3.4 mTorr. Pre-ionization under the coil
(|z| < 0.6 m) is assumed fully ionized at 2 eV.

**Drive:** `drive="coupled"`. The banks (FRC_DESIGN_WINDOW.md) load
self-consistently with the plasma:
- main: eight 6 µF modules, 30 nH each, crowbar 8 nH and about 2 mΩ per module,
  unidirectional start switches;
- bias: NUCTE-like, 2000 µF, 2.4 µH, 8 mΩ.
- Charges are set in vacuum for the stated bias at t = 0 and peak at crowbar.

**Grid:** 24×144, with the flux-conserving field boundary at the coil radius.

**Summary tool:** `python -m examples.design_scan_summary <runs>`.

## Results so far

| Run | Gas | Bias, peak | Coil (R_M) | Charge, bank | Closed FRC from | At 10 µs (r_s, l_s, n̄, Te+Ti, φ) | Later |
|---|---|---|---|---|---|---|---|
| d12_uniform_p020 | D2 | −0.05, 0.20 T | uniform (1.00) | 11.6 kV, 3.2 kJ | 13 µs | still open: reversed to the tube ends | 20 µs: 10.6 cm, 141 cm, 4.7e20, 150 eV, 1.5 mWb |
| d12_uniform_p030 | D2 | −0.05, 0.30 T | uniform (1.00) | 16.2 kV, 6.3 kJ | not by 16 µs | open | open, flux decaying |
| d12_uniform_h2_p020 | H2 | −0.05, 0.20 T | uniform (1.00) | 11.6 kV, 3.2 kJ | 2 µs | 11.4 cm, 86 cm, 6.0e20, 157 eV, 2.0 mWb | 15 µs: 11.1 cm, 110 cm, 5.4e20, 156 eV, 1.85 mWb |
| d12_graded24_p030 | D2 | −0.05, 0.30 T | 0.24 → 0.165 (1.62) | 27.6 kV, 18.3 kJ | 2 µs | 11.4 cm, 58 cm, 9.1e20, 485 eV, 2.75 mWb | 12 µs: l_s 40 cm |

## Findings so far
1. **Without a mirror, D2 forms slowly or not at all.** The reversed field
   stays connected to the tube ends for more than 10 µs. This is NUCTE-III's
   weak-mirror case 1 behaviour, made worse because a uniform coil hugging the
   tube has R_M = 1.00. Two artifact tests (below) rule out the gap floor and
   the tube length as the cause.
2. **H2 forms a closed, elongated FRC on the same 3.2 kJ bank by about 2 µs**
   (r_s ≈ 11 cm, l_s ≈ 0.9–1.2 m, 160 eV). Its lighter ions contract faster,
   and the user's workhorse gas is H2.
3. **A graded coil (R_M 1.6) forms a compact, hot D2 FRC** (l_s 40–60 cm,
   about 480 eV) for about 3× the bank energy.
4. **Cheaper mirrors do not work.** Separately driven end sections reach only
   R_M 1.07 at twice the voltage. Wound mirror-coil pairs need more than 100 kV
   (FRC_DESIGN_WINDOW.md).

## Model-artifact tests (uniform D2 0.3 T)
- **Gap resistivity** (1e-2 Ω·m on near-floor cells): it frees the bias flux
  beyond the coil at first (reversed length 1.31 m instead of 2.17 m at 1 µs),
  but the reversed region regrows to 1.5–1.6 m by 3 µs.
- **Truncated chamber** (gas only to |z| = 0.7 m, exact vacuum beyond, outflow
  ends): the reversed field still reaches the chamber ends.
- **NUCTE-III case 1 with the gap resistivity:** contraction is no earlier
  through 6 µs.

## Caveats
- **2D single-fluid MHD:** no Hall or two-fluid reconnection physics. The end
  reconnection rate in weak mirrors is the least certain piece. In the
  benchmark, case-1 contraction in the model was later than measured (11–15 µs
  vs a contraction onset at 6.5 µs).
- **Pre-ionization is assumed** (frc-ibl.4).
- **Counted dual-energy reset corrections are large at this low fill** (57–72%
  of the coil work by 12–16 µs), but the FRC is insensitive to them. For H2 at
  0.2 T with a 3% floor instead of 1%, from 6 µs on:
  - r_s, n̄ and T agree within 6%, for example 11.53 vs 11.38 cm, 6.4 vs
    6.0e20 m⁻³ and 158 vs 157 eV at 10 µs;
  - l_s agrees within about 15%;
  - only the closure transient before 6 µs differs.
- **Lifetime is not predicted.**
