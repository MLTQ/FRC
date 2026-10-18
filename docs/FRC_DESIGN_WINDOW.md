# First-cut benchtop FRC design window (2026-10-18)

This document gives the fields needed, where and when, to form an FRC by field-reversed
theta pinch in the user's tube. It is a scaling estimate, not a simulation. Its
purpose is to fix the regime that the benchmarked solver (beads frc-ibl.2) will
refine. Expect factor-of-2–3 uncertainty, mostly from the empirical formation
efficiencies listed below.

Reproduce it with `python -m examples.frc_design_window`. The equations are in
[plasma_bench/frc_scaling.md](../plasma_bench/frc_scaling.md).

## Lab envelope used
- Pyrex tube with about 15 cm OD; the bore radius is taken as r_w = 7.1 cm.
- Single-turn sheet main coil hugging the tube: r_c = 8.0 cm, length 0.5 m,
  centred in the ~1 m tube so reconnection has room at both ends.
- Stray (bank, switch, feed) inductance of 30 nH, with 15 and 60 nH as sensitivity cases.
- Gases: H2 as the workhorse, D2 sparingly, N2 for commissioning only.

## Independent checks of the formulas
| Check | Closed form | Independent value |
|---|---|---|
| Rigid-rotor trapped flux, K=1/2/3 (units of π r_s² B_e) | 0.21689 / 0.33125 / 0.38489 | identical to 1e-16 (quadrature) |
| Net flux at the separatrix | 0 | ≤ 4e-17 |
| Classical decay factor τη/(μ0 r_s²), K=1/2/3 | 0.054223 / 0.041406 / 0.032074 | agree to 1e-10 (finite-difference null current) |
| Sheet-coil centre field per ampere | 2.39370e-6 T | 2.39370e-6 (4000 elliptic-integral loops, 3e-9 rel.) |
| Coil inductance | 44.17 nH (Wheeler) | 44.29 nH (flux linkage, 4000 loops; 0.27%) |

These checks verify the arithmetic and geometry. They do not validate the physics assumptions.

## What sets the design (the useful scalings)
- **Temperature** comes from pressure balance: T ∝ B_e² / fill pressure. Doubling
  the field quadruples the temperature at a fixed fill.
- **s̄** (ion gyroradii across the FRC) depends only on the fill and the gas:
  s̄ ∝ √(n/m). It does not depend on the field. Stable formation is reported
  for s̄ of order 1–4; this design targets about 2.
- **Bank energy** depends only on the field and the inductance: E = L I²/2 ∝ B_e².
  A faster rise costs voltage, not energy.
- **The required rise time** comes from the implosion: the main field should
  reach full strength within a few snowplow transit times, t_sp ∝ √ρ_fill / B_e.
  Hydrogen implodes √2 faster than deuterium at the same fill and field. A bank
  fast enough for H2 is therefore also fast enough for D2.
- **Classical lifetime** scales as τ_φ ∝ r_s² T_e^{3/2}. The tube bore is the
  strongest lever: a 10 cm bore doubles τ_φ at 1.6× the energy (see sensitivity).

## Nominal design point
Fill 10 mTorr, B_e = 0.3 T, x_s = r_s/r_c = 0.44.

| | D2 | H2 |
|---|---|---|
| Separatrix radius | 3.55 cm | 3.55 cm |
| ⟨n⟩ (peak) | 2.3e21 (2.5e21) m⁻³ | same |
| T_e + T_i | 88 eV | 88 eV |
| s̄ / r_s/d_i | 1.96 / 5.6 | 2.78 / 7.9 |
| Trapped flux | 0.39 mWb | 0.39 mWb |
| Classical τ_φ (upper bound) | 18 µs | 18 µs |
| Bias field | −62 mT | −62 mT |
| Snowplow time / quarter period | 0.55 / 1.65 µs | 0.39 / 1.17 µs |
| Coil current (single turn) | 125 kA | 125 kA |
| Lossless bank: E, C, V | 583 J, 15 µF, 8.8 kV | 583 J, 7.5 µF, 12.5 kV |
| Minimum PI T_e | 1.5 eV | 1.1 eV |

The single-turn column is the reference case. 125 kA through one switch is
impractical, so the recommended drive splits the current; see the next section.
Build for the H2 rise rate, and the same hardware covers D2 and N2
commissioning.

Stretch option: 0.5 T at 10–20 mTorr D2 gives T_e + T_i of 120–245 eV and a
classical τ_φ of 30–78 µs. It needs 1.6 kJ and 17–25 kV with a single turn.

## Coil drive options (series turns vs parallel modules)
The required rise rate fixes the voltage around each turn of the coil:
V_turn = π r_c² dB/dt at peak. That is 7.4 kV for H2 and 5.3 kV for D2 at the
nominal point, whatever the winding.

- **Series turns:** N turns divide the current by N, but the terminal voltage
  is N × V_turn plus the stray-inductance drop.
- **Parallel modules:** k capacitor/switch modules feed the coil together,
  either at several points on one single-turn coil or through k separate
  single-turn sections. Each module carries I/k at roughly one turn's voltage.

In both cases the stray inductance becomes a smaller share of the total, so the
stored energy falls. The model for synchronously fired modules: each module
sees k N² L_coil + its own 30 nH stray.

Nominal point, 10 mTorr, 0.3 T, lossless:

| Turns × modules | Current per switch | H2 charge V | H2 C per module | D2 charge V | Total energy |
|---|---|---|---|---|---|
| 1 × 1 | 125 kA | 12.5 kV | 7.5 µF | 8.8 kV | 583 J |
| 2 × 1 | 63 kA | 17.4 kV | 2.7 µF | 12.3 kV | 406 J |
| 3 × 1 | 42 kA | 24.0 kV | 1.3 µF | 17.0 kV | 373 J |
| 4 × 1 | 31 kA | 31.0 kV | 0.75 µF | 22.0 kV | 362 J |
| 1 × 2 | 63 kA | 10.0 kV | 4.7 µF | 7.1 kV | 465 J |
| **1 × 4** | **31 kA** | **8.7 kV** | **2.7 µF** | **6.2 kV** | **406 J** |
| 2 × 2 | 31 kA | 16.2 kV | 1.4 µF | 11.4 kV | 376 J |
| 1 × 8 | 16 kA | 8.1 kV | 1.4 µF | 5.7 kV | 376 J |

**Recommendation: four parallel modules on a single-turn coil.** Each module is
2.7 µF at about 9 kV (about 100 J) with its own switch carrying 31 kA. The
total is about 400 J delivered; plan for about 2× stored margin. This divides
the switch current by four without raising the voltage.

Splitting the coil into four 12.5 cm sections, each fed by its own module, also
lets the end sections run harder or carry an extra turn to make the end mirrors.
The cost is that the switches must fire together, with jitter well below the
1.2 µs rise (target ≲ 50 ns). The end sections link less flux, so their
currents will not share evenly; the circuit model (frc-ibl.5) has to resolve
that.

Two series turns on one switch (63 kA, 17 kV) is the alternative if you would
rather have one switch than four. It needs turn-to-turn insulation rated for the
7.4 kV per turn during the µs rise, and the spacing between turns leaves a
small helical ripple in the field near the wall.

Either way of splitting also helps the crowbar. Each crowbar loop has more
inductance for a similar switch and joint resistance: with four modules, each
loop is 207 nH instead of 74 nH. The hold time L/R therefore roughly triples for
the same loop resistance.

## Field timeline (where and when)
| Phase | When (relative to main fire) | Where | What |
|---|---|---|---|
| Vacuum and fill | before shot | tube | Base pressure ≤ 1e-5 Torr (turbo), then a static fill of about 10 mTorr. Glow-discharge clean the Pyrex. |
| Bias | ms before → through reversal | coil length plus ~10 cm beyond each end, ±5% uniform | B_b ≈ −60 mT, opposing the main field. Supply it slowly (or DC) through the main coil via isolation, or from an outer coil if the main coil's eddy time is short. It stores about 25 J. |
| Pre-ionization | tens of µs before main (a scan parameter) | under the coil | Ringing theta PI (HFRC-F used 150 kHz) or a z-discharge, fired with bias at full strength. Target T_e ≳ 2 eV and near-full ionization so the plasma freezes the bias flux. |
| Main rise | 0 → τ_q ≈ 1.2–1.7 µs | main coil | +0.3 T on axis: 125 kA-turns in total, e.g. four 31 kA modules at about 9 kV. The reversal at the wall drives the implosion and the field-line reconnection at the coil ends. |
| Crowbar | at peak current | each module loop | Hold the field. The decay time is L/R. With four modules (207 nH per loop), 50 µs needs R ≤ 4 mΩ per loop and 100 µs needs R ≤ 2 mΩ. A single turn on one module (74 nH) would need 1.5 / 0.74 mΩ. FRC lifetime cannot exceed the field hold. |
| End mirrors | with main | coil ends | Provide for end coils raising the field 1.2–1.5× at the ends to control reconnection and axial contraction. The value comes from the frc-ibl.6 scan, not this estimate. |
| Multipole (later) | through equilibrium | outside the coil | Weak quadrupole windings are reported to suppress the n=2 rotational mode. They only matter once lifetimes exceed tens of µs. |

## Tube bore options (H2, 0.3 T, four parallel modules)
McMaster does not stock a suitable tube. On 2026-10-18 its largest borosilicate
tube was 5.5" OD × 5" ID × 7" long (a pressure-rated sight glass), and its 12"
stock lengths stop at 1" OD. Large bores come from borosilicate process pipe
(DN150 and up) or fused-quartz fabricators. The bore is therefore a free design
choice:

| Tube OD (bore radius) | Fill | r_s | s̄ | T_e+T_i | Classical τ_φ | Module V | Total E |
|---|---|---|---|---|---|---|---|
| 6" (7.1 cm) | 10 mTorr | 3.5 cm | 2.8 | 88 eV | 19 µs | 8.7 kV | 406 J |
| 8" (9.5 cm) | 10 mTorr | 4.8 cm | 3.7 | 87 eV | 33 µs | 10.3 kV | 664 J |
| 8" (9.5 cm) | 5 mTorr | 4.8 cm | 2.6 | 175 eV | 84 µs | 14.6 kV | 664 J |
| 10" (12.0 cm) | 10 mTorr | 6.0 cm | 4.7 | 87 eV | 52 µs | 12.2 kV | 1030 J |
| 12" (14.5 cm) | 10 mTorr | 7.2 cm | 5.7 | 87 eV | 76 µs | 14.2 kV | 1515 J |

s̄ grows with the bore at fixed fill. A larger tube should therefore run at a
lower fill to keep s̄ near 2–3, which also raises the temperature and the
classical lifetime but needs a faster rise. The classical lifetimes remain upper
bounds.

## Sensitivity at the nominal point (D2, 10 mTorr, 0.3 T)
| Varied input | T_e+T_i (eV) | s̄ | classical τ_φ (µs) | bias (mT) | E (J) / V (kV) |
|---|---|---|---|---|---|
| nominal | 88 | 1.96 | 18.5 | −62 | 583 / 8.8 |
| particle retention 0.2 / 0.6 | 176 / 59 | 1.39 / 2.41 | 47 / 11 | −62 | same |
| flux retention 0.2 / 0.6 | 88 | 1.96 | 18.5 | −124 / −41 | same |
| r_s/r_w 0.4 / 0.6 | 59 / 121 | 1.93 / 2.01 | 6.8 / 41 | −40 / −89 | same |
| l_s/l_c 0.4 / 0.8 | 59 / 117 | 2.41 / 1.70 | 11 / 27 | −62 | same |
| T_i/T_e 0.5 / 2 | 88 | 2.41 / 1.70 | 28 / 11 | −62 | same |
| stray L 15 / 60 nH | 88 | 1.96 | 18.5 | −62 | 465/7.1 · 818/12.4 |
| bore 5 / 10 cm | 89 / 87 | 1.38 / 2.78 | 9.3 / 36 | −62 | 412/9.1 · 919/9.5 |

The particle and flux retention are the least certain inputs. They come from
the pre-ionization and reversal physics, which is what the benchmark and the PI
model must pin down.

## Lab implications
1. **The turbo pump is required for H2/D2 FRC work, not optional.** With a
   10 mTorr base pressure and a 10–20 mTorr fill, a third to a half of the gas
   is air. Oxygen and nitrogen line radiation would then dominate the energy
   budget at 10–100 eV. A 1e-5 Torr base keeps base-pressure impurities at 0.1%
   of a 10 mTorr fill.
2. **Crowbar resistance sets the lifetime ceiling.** Splitting the drive into
   modules relaxes it to a few mΩ per loop. A 15 cm tube gives a classical τ_φ
   of about 20 µs at the nominal point, so a crowbar time of 50 µs is sufficient
   for a first device.
3. **Low stray inductance saves both energy and voltage.** Keep each module
   physically close to the coil with low-inductance feeds.
4. **Expect lifetimes well below classical.** Published small FRCs often lose
   flux several times faster than classical resistivity predicts. Treat
   18 µs as an upper bound and plan diagnostics with ≤ 100 ns resolution.
5. **N2 is for commissioning only:** checking bias, PI, timing, the implosion
   signature in excluded flux, and the camera. Nitrogen ionizes to several
   charge states and radiates strongly at these temperatures, so its rows in
   `scan.csv` (which assume Z = 1) are not FRC predictions.

## Experimental anchors
- HFRC-F (Huazhong; [conference paper](https://indico.fusenet.eu/event/28/contributions/101/attachments/115/999/Formation%20studies%20of%20field-reversed%20configurations%20on%20the%20HFRC-F%20device.pdf))
  - Uses the same FRTP sequence as here (bias, ringing PI at 150 kHz, main
    rise, crowbar).
  - Typical shot: bias −15 mT, PI 15 mT, main 35 mT; reported lifetime about
    75 µs at density about 4e20 m⁻³.
  - Reports that PI voltage and timing change the frozen flux, which supports
    treating flux retention as the key empirical input.
  - It is a much larger machine, so it is not a direct size match.
- FAT-CM reports initial FRCs with r ≈ 6 cm and 0.4–0.5 mWb trapped flux at a
  0.4 T compression field, the same scale as the 0.39 mWb nominal point here
  (secondary summary; to be confirmed against the primary paper in frc-ibl.2).

## What this does not establish
This does not establish formation, the reconnection timing, the actual flux or
particle retention, the end dynamics, stability or lifetime. It ignores
radiation, ionization balance and circuit loading. The next step is frc-ibl.2:
choose a documented small FRTP experiment near this regime, reproduce it with
the repository's most complete coil-driven solver path, and only then trust a
design scan of this geometry.
