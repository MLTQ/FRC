# 12-inch design scan: formation runs (frc-ibl.6, 2026-10-19)

## Purpose
Use the benchmarked formation model (NUCTE_BENCHMARK.md) to choose between
build options for the user's 12" device: gas, field, bias and coil shape. Then
state the build specification with its knobs and expected diagnostic signatures.

**What these runs predict:** formation, meaning reversal, end reconnection, axial
contraction, and the size, density, temperature and trapped flux of the
resulting FRC.

**What they do not predict:** lifetime. The flux e-folding times below are 2D
resistive and numerical decay. Real lifetime is set by the n=2 rotational mode
and anomalous transport (frc-ibl.8, FRC_DESIGN_WINDOW.md).

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
  - **uniform** at r = 16.5 cm, or
  - **graded** from r0 at |z| ≤ 0.1 m, linearly down to 16.5 cm at the ends.

**Fill:** D2 at 3.5 mTorr or H2 at 3.4 mTorr. Pre-ionization under the coil
(|z| < 0.6 m) is assumed fully ionized at 2 eV.

**Drive:** `drive="coupled"`. The banks (FRC_DESIGN_WINDOW.md) load
self-consistently with the plasma:
- main: eight 6 µF modules, 30 nH each, crowbar 8 nH and about 2 mΩ per module,
  unidirectional start switches;
- bias: NUCTE-like, 2000 µF, 2.4 µH, 8 mΩ, −0.05 T at main fire.
- Charges are set in vacuum for the bias at t = 0 and the stated peak at
  crowbar.

**Grid:** 24×144, with the flux-conserving field boundary at the coil radius.

**Tools:** `python -m examples.design_scan_summary <runs>` tabulates runs;
`signatures()` in the same module gives the expected diagnostic traces;
`python -m examples.nucte_figure <runs> <png>` plots them.

## Results
"Closed from" is the first time from which the separatrix stays off the wall
and clear of the tube ends for 3 µs. Bias is −0.05 T in every run.

| Run | Gas | Peak | Coil r0 → end (R_M) | Charge, bank | Closed from | At 10 µs: r_s, l_s, n̄ (m⁻³), Te+Ti, φ | At 20 µs | 2D flux e-fold |
|---|---|---|---|---|---|---|---|---|
| d12_uniform_h2_p020 | H2 | 0.20 T | uniform (1.00) | 11.6 kV, 3.2 kJ | 2 µs | 11.4 cm, 86 cm, 6.0e20, 157 eV, 2.02 mWb | 10.9 cm, 83 cm, 5.5e20, 147 eV, 1.76 mWb | 66 µs |
| d12_uniform_h2_p030 | H2 | 0.30 T | uniform (1.00) | 16.2 kV, 6.3 kJ | 6 µs | 9.7 cm, 83 cm, 8.2e20, 250 eV, 1.93 mWb | 9.5 cm, 55 cm, 8.1e20, 222 eV, 1.76 mWb | 130 µs |
| d12_uniform_p020 | D2 | 0.20 T | uniform (1.00) | 11.6 kV, 3.2 kJ | 13 µs | open: reversed to the tube ends | 10.6 cm, 141 cm, 4.7e20, 150 eV, 1.54 mWb | 56 µs |
| d12_uniform_p030 | D2 | 0.30 T | uniform (1.00) | 16.2 kV, 6.3 kJ | not held | open | open; flux e-fold 18 µs | — |
| d12_graded20_p020 | D2 | 0.20 T | 0.20 → 0.165 (1.27) | 15.2 kV, 5.6 kJ | 3 µs | 11.5 cm, 67 cm, 6.5e20, 224 eV, 2.46 mWb | 11.5 cm, 61 cm, 5.7e20, 223 eV, 2.22 mWb | 112 µs |
| d12_graded20_p030 | D2 | 0.30 T | 0.20 → 0.165 (1.27) | 21.3 kV, 10.9 kJ | 2 µs | 9.8 cm, 86 cm, 7.7e20, 325 eV, 2.38 mWb | 10.3 cm, 64 cm, 7.3e20, 336 eV, 2.26 mWb | 183 µs |
| d12_graded24_p020 | D2 | 0.20 T | 0.24 → 0.165 (1.62) | 19.7 kV, 9.3 kJ | 2 µs | 11.1 cm, 98 cm, 5.5e20, 282 eV, 2.84 mWb | 11.8 cm, 67 cm, 5.4e20, 311 eV, 2.58 mWb | 80 µs |
| d12_graded24_p030 | D2 | 0.30 T | 0.24 → 0.165 (1.62) | 27.6 kV, 18.3 kJ | 2 µs | 11.4 cm, 58 cm, 9.1e20, 485 eV, 2.75 mWb | 10.4 cm, 43 cm, 6.8e20, 425 eV, 2.62 mWb | 91 µs |

**Sensitivity runs** (H2, 0.2 T, uniform coil):

| Run | Change from baseline | Closed from | At 10 µs: r_s, l_s, n̄, Te+Ti, φ |
|---|---|---|---|
| d12_uniform_h2_p020_f03 | 3% floor | 6 µs | 11.5 cm, 79 cm, 6.4e20, 158 eV, 2.01 mWb |
| d12_h2_p020_pi1 | pre-ionization at 1 eV | 2 µs | 11.3 cm, 86 cm, 6.1e20, 153 eV, 2.00 mWb |
| d12_h2_p020_pi5 | pre-ionization at 5 eV | 7 µs | 11.5 cm, 89 cm, 5.9e20, 168 eV, 2.03 mWb |

**Model-artifact tests** (uniform coil, D2, 0.3 T, stopped early). Both rule out
the declared gap floor and the tube length as the cause of the slow
uniform-coil D2 formation:
- **Gap resistivity** (1e-2 Ω·m on near-floor cells): it frees the bias flux
  beyond the coil at first (reversed length 1.31 m instead of 2.17 m at 1 µs),
  but the reversed region regrows to 1.5–1.6 m by 3 µs.
- **Truncated chamber** (gas only to |z| = 0.7 m, exact vacuum beyond, outflow
  ends): the reversed field still reaches the chamber ends.

## Findings
1. **H2 forms on the cheapest hardware.** A uniform coil hugging the tube with
   H2 forms a closed, elongated FRC: 11 cm × 0.8–1.1 m, about 150 eV, 2 mWb,
   by 2–6 µs on the 3.2 kJ first-light bank.
   - At 0.3 T (6.3 kJ) it is hotter (220–250 eV) and denser, but closes later
     (6 µs).
2. **D2 needs a mirror.** A uniform coil hugging the tube has R_M = 1.00. With
   D2, the reversed field stays tied to the tube ends: it closes only at 13 µs
   at 0.2 T and never holds at 0.3 T. NUCTE-III's weak-mirror case 1 behaves
   the same way.
   - **Mild graded coil** (centre radius 20 cm, R_M 1.27): D2 forms by 2–3 µs.
     At 0.2 T (5.6 kJ) it gives 11.5 cm × 0.6–0.7 m, 220 eV, 2.2–2.5 mWb. At
     0.3 T (10.9 kJ), 330 eV.
   - **Stronger graded coil** (24 cm, R_M 1.62): more flux (2.6–2.8 mWb) and
     hotter, for 1.7× more energy again.
   - **Other mirrors:** driving short end sections harder (R_M 1.07 at 2×
     voltage) and wound mirror coils (> 100 kV) are not practical
     (FRC_DESIGN_WINDOW.md).
3. **Robustness.**
   - The 1% vs 3% gap floor agrees within about 6% after closure, and l_s
     within about 15%.
   - Pre-ionization temperature from 1 to 5 eV changes the FRC by ≤ 10%.
   - Trapped flux is set by the bias, not the pre-ionization temperature, for a
     fully ionized pre-ionization. Partial ionization is not modeled
     (frc-ibl.4).
4. **More mirror holds flux longer in 2D:** the mild graded coil at 0.3 T e-folds
   in 183 µs. This is a 2D trend only; it does not predict lifetime.

## Recommended build (formation)
Common to every option:
- 12" ID tube, 2.2 m long.
- Single-turn coil, 1.2 m long, all straps on one collector.
- Eight 6 µF modules rated at least 20 kV (the bank is sized for the upgrade
  path), with unidirectional start switches.
- One crowbar per module, ≤ 8–20 nH and about 2 mΩ each.
- One trigger fanned out over matched cables (±100–200 ns is fine).
- Bias bank about 2000 µF, about 2.4 µH isolation, about −2.3 kV for −0.05 T,
  fired about 100 µs before main.
- Pre-ionization bank about 10 µF at about 3.5 kV, ringing at 150 kHz, fired
  tens of µs before main.

| Option | Coil | Gas, fill | Main charge | Expected FRC | Notes |
|---|---|---|---|---|---|
| **A, first light (recommended start)** | uniform, r = 16.5 cm | H2, 3.4 mTorr | 11.6 kV (3.2 kJ) for 0.2 T | 11 cm × 0.8–1.1 m, 150 eV, n̄ 6e20, 2.0 mWb | Cheapest coil and bank. H2 is the workhorse gas. Raising the charge to 16.2 kV gives 0.3 T and about 230 eV. |
| B, D2 or hotter | mild graded, 20 cm → 16.5 cm (taper from \|z\| = 0.1 m) | D2, 3.5 mTorr (H2 also) | 15.2 kV (5.6 kJ) for 0.2 T; 21.3 kV (10.9 kJ) for 0.3 T | 11.5 cm × 0.6–0.7 m, 220 eV (0.2 T); 10 cm × 0.6–0.9 m, 330 eV (0.3 T) | About 1.7× the bank energy. Needed for D2. More trapped flux (2.2–2.5 mWb). |

The coil shape is the one decision that is hard to change later:
- A uniform coil limits D2 operation.
- A mild graded coil (20 cm centre) costs 1.7× the bank energy but works for both
  gases.
- Option A's bank (8 × 6 µF rated ≥ 20 kV) drives option B at 15.2 kV.

## Expected diagnostic signatures (midplane, 8–20 µs)
From `signatures()`. The interferometer is side-on through the axis; the wall
probes are just inside the tube.

| Option | Wall B_z (B-dot) | Excluded-flux radius (flux loop) | ∫n dl (interferometer) | FRC length (camera) | Peak coil current | Capacitor voltage held after crowbar |
|---|---|---|---|---|---|---|
| A: H2, uniform, 0.2 T | 0.22–0.26 T | 10.9–11.8 cm | 1.3–1.6e20 m⁻² | 0.64–1.19 m | 228 kA at 3 µs | −2.8 kV |
| A at 0.3 T | 0.31–0.36 T | 9.5–10.2 cm | 1.6–2.0e20 m⁻² | 0.55–0.86 m | 315 kA at 5 µs | −4.0 kV |
| B: D2, graded, 0.2 T | 0.26–0.32 T | 11.3–12.4 cm | 0.95–2.3e20 m⁻² | 0.52–1.01 m | 255 kA at 4 µs | −3.3 kV |
| B at 0.3 T | 0.35–0.42 T | 10.1–11.2 cm | 1.2–2.2e20 m⁻² | 0.40–0.86 m | 394 kA at 3 µs | −4.7 kV |

**Formation signature:**
- The excluded-flux radius drops below the wall within 1–2 µs, then stays near
  10–12 cm.
- The camera shows the luminous column shortening from about the coil length to
  under 1 m between 3 and 8 µs, then bouncing axially once or twice.
- A uniform-coil D2 shot that fails shows a column that never shortens.

## Caveats
- **2D single-fluid MHD:** no Hall or two-fluid reconnection physics. The end
  reconnection rate in weak mirrors is the least certain piece. In the
  benchmark, the model's weak-mirror contraction was later than measured.
- **Pre-ionization is assumed fully ionized,** so the ionization fraction is
  untested (frc-ibl.4).
- **Counted dual-energy reset corrections are large at this low fill** (57–72%
  of the coil work by 12–16 µs). They live in the gap, and the FRC is
  insensitive to them (floor test above).
- **Lifetime is not predicted:** see the n=2 section of FRC_DESIGN_WINDOW.md.
  Start without multipoles, measure the onset, then add the quadrupole.
- **Resolution:** 24×144. The benchmark showed formation converged between
  16×128 and 32×256 for NUCTE. The design geometry has not had its own
  convergence run.
