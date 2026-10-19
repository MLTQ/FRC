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

## 12-inch bore for maximum lifetime (user target, 2026-10-18)
**Geometry:**
- Bore: 12" ID, r_w = 15.2 cm.
- Coil: single turn, r_c = 16.5 cm, 1.2 m long, fed by eight parallel modules.
  A 0.5 m coil would make a squat FRC; 1.2 m keeps l_s/(2r_s) ≈ 4.7, inside
  the elongated-FRC relations used here.
- Tube: about 2–2.4 m long.

The fill is chosen to hold s̄ at 2.5 or 3.5 (`large_bore.csv`).

| Gas, fill | B_e | s̄ | T_e+T_i | Classical τ_φ | Loop V per turn | Module V | Total E |
|---|---|---|---|---|---|---|---|
| D2, 3.5 mTorr | 0.2 T | 2.5 | 111 eV | 113 µs | 7.9 kV | 8.3 kV | 1.6 kJ |
| H2, 3.4 mTorr | 0.2 T | 3.5 | 113 eV | 116 µs | 11.3 kV | 11.9 kV | 1.6 kJ |
| D2, 6.8 mTorr | 0.3 T | 3.5 | 128 eV | 141 µs | 12.8 kV | 13.4 kV | 3.7 kJ |
| D2, 3.5 mTorr | 0.3 T | 2.5 | 250 eV | 354 µs | 17.9 kV | 18.7 kV | 3.7 kJ |
| H2, 3.4 mTorr | 0.3 T | 3.5 | 255 eV | 363 µs | 25.5 kV | 26.7 kV | 3.7 kJ |
| D2, 3.5 mTorr | 0.4 T | 2.5 | 445 eV | 799 µs | 31.7 kV | 33.2 kV | 6.6 kJ |
| H2, 1.7 mTorr | 0.5 T | 2.5 | 1390 eV | 3.9 ms | 99 kV | 104 kV | 10.3 kJ |

**What limits the lifetime:**
1. **Classical flux decay.** At fixed s̄, τ_φ ∝ B_e³ r_s⁵ / (m_i^{3/2} s̄³).
   Field and bore dominate; at equal fill and field, H2 and D2 reach the same
   temperature and the same classical τ_φ.
2. **Formation voltage.** The loop voltage per turn grows about as B_e² at
   fixed fill, because the field must also rise faster. This is why
   field-reversed theta-pinch formation gets hard in large, high-s̄ devices
   (LSX needed "elaborate formation controls"). D2 implodes more slowly than H2
   at the same fill, so it needs about 70% of the H2 voltage. That makes D2
   the better large-bore gas.
3. **Measured confinement is below classical.** Reported results:
   - LSX: up to 0.5 ms for s̄ ≤ 4.
   - FRX-C: over 300 µs with weak quadrupole fields suppressing the n=2
     rotational mode, with particle confinement scaling as R².

   Both machines are larger than this bore. Scaling them down by R² suggests
   roughly **50–150 µs** is realistic here. The table's classical values are
   ceilings.

   A closer anchor is NUCTE-III (NUCTE_BENCHMARK.md): coil radius 0.17 m, D2
   at 10 mTorr. It measured 50–65 µs, which is about 0.5–0.9 of this model's
   classical estimate for its geometry.
4. **The n=2 rotational instability** ends theta-pinch FRCs unless suppressed.
   Quadrupole windings are required for this target (beads frc-ibl.8).
5. **Field hold.** Hold at least 3× the target lifetime: 300 µs with
   L_mod = 671 nH needs R ≤ 2.2 mΩ per module loop.

**Recommended 12" baseline:**
- **Gas and fill:** D2 at about 3.5 mTorr, s̄ ≈ 2.5.
- **Modules:** eight, each about 6 µF and rated to at least 20 kV. The circuit
  model below supersedes the earlier 15 kV figure.
- **First light:** charge to 11.6 kV for 0.2 T, starting from a −0.05 T bias
  (3.2 kJ; circuit model below). This meets the rise criterion; classical
  τ_φ ≈ 113 µs.
- **Upgrade:** charge the same capacitors to 16.2 kV for 0.3 T (6.3 kJ).
  - The lossless single-turn estimate gave 8.3 kV and 12.5 kV. It omitted three
    things: the field must swing up from the bias, the bias circuit diverts some
    current, and the equal-flux coil is 77 nH.
  - The quarter period stays at 3.1 µs, which is 4.5 snowplow times instead of
    3. Whether that slower relative rise still forms well is for frc-ibl.6 to
    determine.
- **Gas use:** D2 consumption is negligible. A static fill of a
  12" × 2.2 m tube at 7 mTorr is about 1.4 mL at STP per shot.
- **Vacuum:** base pressure of about 1e-6 Torr, which keeps base-gas
  impurities near 0.1% of a 3.5 mTorr fill. A turbo and a clean or baked tube
  are required.
- **Radiation:** D2 at a few hundred eV will produce some D-D neutrons; this is
  for the radiation consultant.

## Bank, crowbar, timing and mirror for the 12" baseline (circuit model, frc-ibl.5)
Computed with `plasma_bench/coil_circuit.py` (verified in its companion doc) by
`python -m examples.frc_bank_design`.

**Assumed hardware:**
- **Coil:** single turn, r = 16.5 cm, 1.2 m long, 24 straps of 5 cm on one
  collector.
- **Main bank:** eight modules of 6 µF with 30 nH stray each, fired together.
- **Bias bank:** NUCTE-like, 2000 µF, 2.4 µH, 8 mΩ, fired so it peaks at main
  fire.
- **Charge voltages** are set so the vacuum field at the coil centre is the bias
  value at t = 0 and the stated peak at crowbar.

**Coil:**
- Equal-flux inductance 77.3 nH; centre field 0.903 mT per kA.
- A uniform single-turn coil has no mirror: its field is maximal at the
  midplane.

| Peak (from bias) | Main charge | Bank energy | Bias charge | Peak per switch | dI/dt per switch | Crowbar |
|---|---|---|---|---|---|---|
| 0.2 T (−0.05 T) | 11.6 kV | 3.2 kJ | −2.3 kV | 36 kA | 18 kA/µs | 3.05 µs |
| 0.3 T (−0.05 T) | 16.2 kV | 6.3 kJ | −2.3 kV | 50 kA | 26 kA/µs | 3.05 µs |
| 0.3 T (−0.07 T) | 17.1 kV | 7.0 kJ | −3.2 kV | 53 kA | 27 kA/µs | 3.05 µs |

### Crowbar and start switches
**Use unidirectional start switches** (ignitrons, or spark gaps with a series
diode stack). Without them, the 6 µF capacitors keep ringing through the crowbar
after it fires and pull the held field down. The circuit shows the same effect
on NUCTE's bank: −20% within 6 µs.

**Keep the crowbar loop inductance small next to the 77 nH coil.** When the
crowbar fires, the coil keeps only L_coil/(L_coil + L_crowbar) of its current.

| Crowbar L per module (8 in parallel) | Field held after 2 µs | After 10 µs | After 25 µs | Capacitor reversal |
|---|---|---|---|---|
| 8 nH | 98% | 95% | 89% | 25% of charge |
| 20 nH | 96% | 93% | 87% | 28% |
| 40 nH | 92% | 89% | 84% | 33% |

- **Hold time:** a 300 µs L/R hold needs about 2.1 mΩ per module crowbar loop,
  switch and joints included.
- **Capacitor rating:** rate the capacitors for about 30% voltage reversal.

### Timing
**Common collector (all eight modules on one coil):**
- The coil's equal-flux current distribution fixes the field shape, so skew only
  reshapes the total current.
- One module firing 50, 200 or 500 ns late changes the centre field by 0.16%,
  0.6% or 1.8% during the rise.
- It moves that switch's peak current by up to about ±7%.

**Separate axial sections (one module each):** skew becomes axial asymmetry.
With the four modules on z > 0 late:

| Skew | Asymmetry at 0.5 µs | Asymmetry at 1 µs |
|---|---|---|
| 10 ns | 3.6% | 1.1% |
| 50 ns | 18% | 5.5% |
| 200 ns | 72% | 23% |

**What the asymmetry does to an FRC** (NUCTE-III case 6, two halves fired 100 or
300 ns apart, formation runs):
- trapped-flux centroid displacement ≤ 0.5 or 1.5 cm;
- left/right flux asymmetry ≤ 2.4% or 6%;
- both small against the 40–60 cm FRC.

**Recommendation:** a common collector with a single trigger fanned out to the
eight switches over matched cables. ±100–200 ns of jitter is tolerable, so no
phase-locked timing is needed.

### End mirror
NUCTE-III makes its mirror by stepping the coil radius down at the ends. Its
case 6 (R_M 1.8, field index 0.3) contracts by 6–7 µs. Its case 1 (R_M 1.1)
contracts only by 15–25 µs, yet lived longer (65 vs 50 µs, n=2 onset 28 vs
22 µs).

**Option 1: graded single-turn coil.** The radius widens from 16.5 cm at the
ends to r0 at |z| ≤ 0.1 m; values for 0.3 T at the centre on the same 8 × 6 µF
bank:

| r0 | R_M | Field index | Coil | Charge | Energy |
|---|---|---|---|---|---|
| 0.20 m | 1.27 | 0.13 | 93 nH | 21 kV | 11 kJ |
| 0.22 m | 1.44 | 0.20 | 102 nH | 24 kV | 14 kJ |
| 0.24 m | 1.61 | 0.27 | 110 nH | 28 kV | 18 kJ |

**Option 2: separate mirror coil pairs.** N-turn coils just beyond the coil
ends, as one series circuit. For R_M 1.3–1.5 they need 10–15% of the main
current per turn in 2–4 turns. That is a 4–13 µH circuit, needing more than
100 kV for a µs rise, so it is impractical.

**Decision:** whether the mirror is worth up to 3× the bank energy is a
formation question for frc-ibl.6. NUCTE's data suggest a weak mirror forms more
slowly but is not worse for lifetime.

### Pre-ionization bank (ringing θ-PI)
**Model.** The PI bank rings on the same 77 nH coil, with the bias bank's
2.4 µH branch in parallel. The wall electric field at the first current peak is
E_θ = ω B_PI r_w/2.

| Frequency | PI stray | C | 20 V/cm at the wall | 40 V/cm at the wall |
|---|---|---|---|---|
| 150 kHz | 40 nH | 9.8 µF | B_PI 0.028 T: 3.5 kV, 60 J, 31 kA | B_PI 0.056 T: 7.0 kV, 240 J, 62 kA |
| 150 kHz | 80 nH | 7.3 µF | 4.7 kV, 80 J | 9.4 kV, 320 J |
| 100 kHz | 40 nH | 22 µF | B_PI 0.042 T: 3.5 kV, 140 J, 46 kA | B_PI 0.084 T: 7.0 kV, 550 J, 92 kA |

For comparison, NUCTE-III used 11.25 µF ringing at 124 kHz, up to 135 kA.

**Constraints:**
- **Keep the PI field below the bias.** The ringing field adds to the bias, so
  |B_PI| < |B_bias| (0.05 T) keeps the field from reversing during PI. That
  limits a −0.05 T bias to roughly 20–30 V/cm at 150 kHz. A stronger PI needs a
  stronger bias or a higher frequency.
- **Conductivity is not the hard part.** The minimum PI electron temperature for
  the plasma column to hold the bias flux through the main rise (skin time ≥ 10
  quarter periods) is 0.56 eV here.
- **What the formation runs assume:** full ionization at 2 eV. Partial
  ionization is not modeled (frc-ibl.4).

### n=2 rotational stabilization (literature, frc-ibl.8; not simulated)
Theta-pinch FRCs spin up, and the n=2 rotational instability ends them. NUCTE-III
saw its onset at 22–28 µs, and its lifetimes (50–65 µs) are set by it. The
axisymmetric solver cannot model it. These numbers come from experiments:

- **Threshold:** a quadrupole field of about 7% of the external confinement field
  suppressed n=2 completely in FRX-C/T with straight coils, and 9% with helical
  coils (Tuszewski et al., Phys. Fluids 29, 852 (1986)).
- **First demonstration:** suppression with a quadrupole field "much smaller than
  the axial confinement field at the separatrix", about 2.5× below the MHD
  prediction (Ohi et al., PRL 51, 1042 (1983)).
- **Payoff:** FRX-C lifetimes exceeded 300 µs once n=2 was stabilized.
- **Penalty:** quadrupole-stabilized FRCs shrink axially and lose particles
  faster. Hexapole and octupole fields are reported with much smaller penalties
  (JPSJ 63, 2845 (1994); Fusion Technol. 27 (1995)), but they fall off faster
  with radius.

**Straight quadrupole for the 12" baseline.** Four axial bars at a = 0.20 m,
just outside the main coil, with currents ±I alternating:
- Field inside: |B| = (2μ0 I/(π a))(r/a). This leading-order form is within
  1–7% of the exact line-current sum out to the wall.
- For 7% of B_e = 0.35 T, each bar carries:
  - about 15 kA if the 7% is set at r_s = 8 cm;
  - about 8 kA if it is set at the wall.

  The papers quote the field "at the separatrix"; budget for the larger number.
- Use several turns per bar to bring the current down.

**Timing.** Establish the quadrupole slowly, over a ms or more, before the main
bank fires:
- A steady m=2 field induces no currents in the single-turn sheet coil, which
  then keeps it frozen through the µs compression.
- A quadrupole pulsed on µs scales would be screened by the sheet coil unless
  it is slotted.

**Recommendation:**
- Start without multipoles to measure the n=2 onset and rotation (camera, end-on
  view).
- Add the straight quadrupole once lifetimes approach the onset time.
- Consider an octupole if particle loss with the quadrupole is a problem.

### Plasma loading of the bank
**For the 12" baseline:** the bank's stray inductance (30 nH / 8 = 3.75 nH) is
5% of the coil's, so the plasma's flux exclusion barely loads it. The early
field drops by ≲ 2–3%.

**For NUCTE-III** (27 nH stray against a 58–69 nH coil), the effect is large:
- the coupled formation run holds the coil's vacuum current 22% below the
  prescribed waveform at 1 µs;
- the gap is 5–8% at 2–3 µs.

Keep strays small. That is what lets the circuit, rather than the plasma, set
the field.

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
