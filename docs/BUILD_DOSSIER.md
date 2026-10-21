# Build dossier: 12-inch benchtop FRC (2026-10-20)

This page collects the electrical and timing requirements for building the
device: what each coil must carry, when, and what bank produces that. It does
not cover hardware engineering, component selection or safety, which are yours.

**Where the numbers come from:**
- the verified circuit model (`plasma_bench/coil_circuit.py`), with plasma
  loading;
- the benchmarked formation runs (NUCTE_BENCHMARK.md, DESIGN_SCAN.md,
  EXPERIMENT_PLAYGROUND.md);
- for the quadrupole only, the literature.

Lines marked *pending* wait on simulations now running.

## 1. Configuration at a glance
| Item | Option A (first light, recommended) | Option B (D2 / hotter) |
|---|---|---|
| Tube | 12" ID borosilicate, 2.2 m | same |
| Main coil | single turn, 1.2 m long, 24 × 5 cm straps on one collector, uniform r = 16.5 cm | same straps, graded: r = 20 cm for \|z\| ≤ 0.1 m, tapering to 16.5 cm at the ends (mirror 1.27) |
| Coil inductance (equal flux) | 77 nH | 93 nH |
| Field per kA (centre) | 0.90 mT/kA | — |
| Gas, fill | H2, 3.4 mTorr (0.45 Pa) | D2, 3.5 mTorr (H2 also works) |
| Field sequence | bias −0.05 T → main peak 0.2 T (0.3 T is hotter but shorter-lived in 2D, §2) | same |
| Expected FRC (2D model) | 11 cm × 0.8–1.1 m, Te+Ti ≈ 150 eV (0.2 T), ≈ 230 eV (0.3 T), 2.0 mWb | 11.5 cm × 0.6–0.7 m, 220 eV (0.2 T), 330 eV (0.3 T), 2.2–2.5 mWb |

**Annular variant** (second tank, for flux-core and toroidal-field
experiments): the same tube and coil, plus a coaxial 3" OD insulating inner tube.
It costs about 40% of the trapped flux at formation (EXPERIMENT_PLAYGROUND.md).

## 2. Main bank (theta coil)
Eight modules in parallel on the common collector, fired together:
- each module 6 µF, about 30 nH stray (3.75 nH for all eight);
- one start switch and one crowbar per module.

| Case | Charge | Bank energy | Peak per switch | dI/dt per switch | Total coil current (vacuum) | Rise (crowbar) | Capacitor reversal |
|---|---|---|---|---|---|---|---|
| A, 0.2 T | 11.6 kV | 3.2 kJ | 36 kA | 18 kA/µs | 221 kA | 3.05 µs | −2.9 kV (25%) |
| A, 0.3 T | 16.2 kV | 6.3 kJ | 50 kA | 26 kA/µs | 332 kA | 3.05 µs | −4.0 kV |
| B, 0.2 T | 15.2 kV | 5.6 kJ | 43 kA | 20 kA/µs | 266 kA | 3.33 µs | −3.5 kV |
| B, 0.3 T | 21.3 kV | 10.9 kJ | 60 kA | 29 kA/µs | 399 kA | 3.33 µs | −4.8 kV |

**0.3 T does not buy lifetime in 2D.**
- Run to 40 µs, option A at 0.3 T matches 0.2 T up to 20 µs (1.76 mWb).
- It then collapses: 0.24 mWb at 40 µs, against 1.11 mWb at 0.2 T.
- The hotter plasma loses 2.2 × the energy along open field lines
  (EXPERIMENT_PLAYGROUND.md, Lifetime levers).
- The earlier "130 µs e-fold" was fitted over 8–20 µs only.
- Use 0.3 T for hotter formation studies, not for lifetime.

**Ratings implied by the table:**
- Capacitors: at least 20 kV for option A with the 0.3 T charge. Plan for 25 kV if
  option B at 0.3 T is on the roadmap. Rate for about 30% voltage reversal.
- Start switches: unidirectional (ignitron, or a spark gap with a series diode
  stack). Without them the capacitors keep ringing through the crowbar and pull the
  held field down; on NUCTE's bank the circuit model shows −20% within 6 µs.
- Crowbars: one per module, ≤ 8–20 nH and about 2.1 mΩ per module loop
  (switch and joints included) for a 300 µs L/R hold.
  - One-way (conducting only freewheeling current) is preferred. Two-way
    crowbars short any later bank on the same coil, such as a re-boost.
  - The held field after crowbar: 98% / 95% / 89% at 2 / 10 / 25 µs with 8 nH
    crowbars; 96% / 93% / 87% with 20 nH.
- Triggering: one trigger fanned out over matched cables. On a common collector
  ±100–200 ns of jitter changes the field by under 1%, so no phase-locked timing
  is needed.
- Magnetic pressure on the coil: 16 kPa (2.3 psi) at 0.2 T, 36 kPa (5.2 psi)
  at 0.3 T.

## 3. Bias bank (reverse field, before main)
The bias field is −0.05 T at main fire, opposing the main field, through the
same coil via an isolation inductor. The isolation keeps the main pulse out of
the bias bank.

| Version | Bias bank | Isolation L | Charge | Energy | Peak current | Fire before main |
|---|---|---|---|---|---|---|
| Simulated baseline (NUCTE-like) | 2000 µF | 2.4 µH, 8 mΩ | −2.3 kV (−2.8 kV for option B) | 5.3 kJ (7.7 kJ for B) | 64 kA | 103 µs (bias peaks at main fire) |
| Cheaper alternative (circuit-checked) | 3.8 mF | 1.0 µH | −1.3 kV | 3.2 kJ | 76 kA | 88 µs |

With the cheaper alternative, the lower isolation diverts more main current
into the bias branch, so the main charge rises from 11.6 to 11.8 kV for the same
0.2 T. Most of the bias energy sits in the isolation inductor.

A slow, multi-turn bias coil outside the main coil is another route. It only
works if the single-turn coil's eddy time is short (a slotted sheet); not
modelled.

## 4. Pre-ionization (ringing theta PI)
**Bank:** about 10 µF at 3.5 kV, with about 40 nH stray.
- It rings at about 150 kHz on the main coil, in parallel with the bias branch.
- 60 J, 31 kA peak.
- It gives B_PI ≈ 0.028 T and about 20 V/cm at the wall.

**Timing:** fire with the bias near its peak, tens of µs before main. Treat the
delay (about 10–50 µs) as a scan parameter.

**Constraints:**
- Keep |B_PI| below |B_bias| (0.05 T) so the PI ringing never reverses the
  field. That caps the PI at about 20–30 V/cm at 150 kHz with this bias. A
  stronger PI needs a stronger bias or a higher frequency (100 kHz needs
  22 µF).
- A plasma at Te ≥ 0.6 eV already holds the bias flux through the rise, so
  conductivity is not the limit.

**What the model assumes:** the formation runs assume full ionization at 2 eV.
They are insensitive to the temperature (1–5 eV changes the FRC by ≤ 10%).

**Not modelled:**
- the ionization fraction (frc-ibl.4);
- the PI voltage and timing, which experiments report change the trapped flux.

Expect to tune these on the bench. N2 is useful for commissioning bias, PI and
timing, but its plasma radiates strongly and is not an FRC prediction.

## 5. Timeline of one shot (t = 0 at main fire)
| When | What |
|---|---|
| ms–s before | Turbo to ≤ 1e-6 Torr base, glow-discharge clean, static fill of 3.4 mTorr H2 |
| ms before (upgrade) | Quadrupole on, slowly (§6) |
| −103 µs | Bias bank fires (peaks at t = 0) |
| about −50 to −10 µs | PI bank fires (scan) |
| 0 | Main bank fires: 8 switches, one trigger |
| about 1–2 µs | Radial implosion; excluded-flux radius drops below the wall |
| 3.05 µs | Crowbars close at capacitor voltage zero (field held, 300 µs L/R) |
| 2–6 µs | Closed FRC forms (end reconnection, axial contraction) |
| 3–8 µs | Column shortens to under 1 m, then bounces axially once or twice |
| 4 µs (annular only) | Flux-core bank fires (§7) |
| 15 µs (optional) | Re-boost bank (§7) |
| tens of µs | n=2 rotation onset expected unless the quadrupole is on |
| 50–150 µs | Realistic FRC lifetime range (literature scaling, with n=2 suppressed) |

Record diagnostics over at least 0–200 µs with ≤ 100 ns resolution.

## 6. Stability: quadrupole for the n=2 rotational mode
Literature-based (Tuszewski 1986; Ohi 1983); not simulated in 2D.
- **Field:** a straight quadrupole of about 7% of the external field (9% for
  helical windings). FRX-C exceeded 300 µs once n=2 was stabilized.
- **Bars:** four axial bars at r = 0.20 m, just outside the main coil, with
  alternating currents.
  - For 7% of 0.35 T, about 15 kA per bar if set at the separatrix (r ≈ 8 cm),
    or about 8 kA if set at the wall.
  - Budget the 15 kA. Use several turns per bar, e.g. 10 turns × 1.5 kA.
- **Timing:** establish it slowly, over ms, before the main bank. A steady m=2
  field induces no current in the single-turn sheet coil, which then holds it
  through the µs compression. A µs-pulsed quadrupole would be screened.
- **Particle cost (test-particle orbits, EXPERIMENT_PLAYGROUND.md "Particle
  orbits"):**
  - Extra ion loss at 0.2 T over 10–30 µs, against 5% without the
    quadrupole: about 24 points at 7% B_e set at the wall (5.5 kA), and
    33 points set at the separatrix (7.5 kA).
  - The lost ions strike the tube wall at |z| ≈ 0.5–1.1 m, in two sectors
    180° apart in each half of the tube. The pattern is rotated 90° between
    the two halves.
  - These are vacuum-field upper bounds.
- **Plan:**
  1. First shots without it, to measure the n=2 onset (end-on camera).
  2. Then add it.
  3. Consider an octupole if quadrupole particle loss is large.

## 7. Optional coils and banks (playground results)
| Coil | Geometry | Bank (as simulated) | Peak current | What it does |
|---|---|---|---|---|
| Re-boost (main coil) | the main coil itself | 8 × 6 µF at 6 kV (0.9 kJ), fired at 15 µs | 38 kA per module | Compresses the FRC (smaller, denser, hotter); adds no trapped flux. A smaller step (8 × 6 µF at 3 kV, 216 J, 15 µs, no crowbar of its own) holds the external field of the 2 kV flux-core case within 2% of the control, but leaves the FRC unchanged. **Not needed for lifetime.** Needs one-way main crowbars. |
| **Flux core, recommended** (annular tank) | 40 turns, r = 2.5 cm, 0.8 m long, inside the inner tube; 4.8 µH plus ~1 µH lead; M = 71 nH to the main coil | **300 µF at 1.0 kV (150 J)**, fired at **4 µs**, main-field polarity, crowbar at its voltage zero | **5.9 kA** circuit peak, ~65 µs after firing (≈ 4 kA in the coil at 40 µs) | Holds the trapped flux flat at 1.33–1.41 mWb from 10 to 40 µs (control 1.24 → 0.71). The FRC still shortens to E ≈ 1.2 by 40 µs, a tilt risk. 0.7–1.0 kV brackets flat flux. The start switch must stay on through a brief current reversal at turn-on: a switch that opens at the first current zero lost the 0.7 kV core at 4.6 µs, and the 1 kV core came within 23 A of opening. |
| Flux core, 2 kV | as above | 300 µF at 2 kV (0.6 kJ), fired at 4 µs | 11.4 kA circuit peak | Overdrives: φ rises to 2.8 × the control by 40 µs, and the FRC turns short and fat (E 1.0). Holding the external field with a re-boost does not change that: the shortening comes from energy loss at constant flux, not from field sag. |
| Slow mirror (bias time) | 4-turn pairs at z = ±0.635 m, r = 17 cm, in series; ≈ 15.5 µH | 200 µF at 3 kV (0.9 kJ), 0.5 µH lead, fired at −100 µs with the bias, crowbarred at its peak | 11.7 kA, falling to 4.9 kA through the main rise | **Not recommended.** The parallel main straps spread its flux over the bore, cancelling the bias unless the bias bank rises to −4.4 kV (19 kJ instead of 5.3 kJ). At peak, B on axis at the coil end rises only from 0.72 to 0.88 of the centre field. The FRC bounces less but has no more flux by 30 µs. |
| Toroidal field (annular tank) | axial rod inside the inner tube, coaxial return cage outside the main coil; ≈ 1.3 µH | estimate: about 3 mF at about 0.4–0.5 kV (≈ 0.3 kJ) for a ~100 µs rise | 20 kA → B_θ ≈ 0.05 T at r = 8 cm; 40 kA → 0.1 T | Twists the field lines (magnetic shear); too much makes it a spherical torus. *Pending:* toroidal-field session. |
| End plugs | 1-turn pairs at z = ±0.635 m, r = 17 cm | 20 µF at 15 kV (2.25 kJ), fired at 3 µs | 69 kA | **Not recommended.** Fired after the crowbar they only add field beyond the coil ends (0.14 → 0.25 T at z = 0.6 m), and the flux-conserving main coil drops its current from 150 to 86 kA. Trapped flux is unchanged. A mirror has to come from the coil shape (option B); slow coils energized with the bias do not form one either (row above). |
| Core punch | 10 turns, r = 2.5 cm, 0.3 m, inside the inner tube | 1.6 µF at 19 kV (0.29 kJ), 2 µs rise, fired at 8 µs | 24.7 kA (1 T inside) | Pushes the FRC's inner edge from 4.5 to 6.8 cm, then splits the FRC. A disruption experiment, not a lifetime tool. |

## 8. Energy per shot (option A, 0.2 T)
| Bank | Energy |
|---|---|
| Main | 3.2 kJ |
| Bias, baseline | 5.3 kJ (3.2 kJ with the alternative) |
| PI | 0.06 kJ |
| **Total stored** | **≈ 8.6 kJ** (≈ 6.5 kJ with the alternative bias) |

Optional banks add:
- re-boost 0.9 kJ;
- flux core 0.15 kJ (1 kV, recommended; 0.6 kJ at 2 kV);
- toroidal rod ≈ 0.3 kJ;
- quadrupole: small; slow ms supply.

At 0.3 T the main bank becomes 6.3 kJ.

## 9. Expected diagnostic signals (midplane, 8–20 µs; DESIGN_SCAN.md)
| Case | Wall B_z (B-dot) | Excluded-flux radius | ∫n dl (side-on interferometer) | FRC length (camera) | Peak coil current |
|---|---|---|---|---|---|
| A, 0.2 T | 0.22–0.26 T | 10.9–11.8 cm | 1.3–1.6e20 m⁻² | 0.64–1.19 m | 228 kA at 3 µs |
| A, 0.3 T | 0.31–0.36 T | 9.5–10.2 cm | 1.6–2.0e20 m⁻² | 0.55–0.86 m | 315 kA at 5 µs |
| B, 0.2 T | 0.26–0.32 T | 11.3–12.4 cm | 0.95–2.3e20 m⁻² | 0.52–1.01 m | 255 kA at 4 µs |

**A failed shot** (for example D2 on the uniform coil) shows a column that
never shortens.

**Hall physics** should leave a small out-of-plane (B_θ) field near the coil
ends in a four-lobed pattern, which a B_θ-oriented probe could see. *Pending:*
the Hall session will give its size.

## 10. Vacuum, gas and radiation
- **Vacuum:** a turbo is required. Aim for ≤ 1e-6 Torr base, which keeps
  impurities near 0.1% of a 3.4 mTorr fill. Glow-discharge clean (H2) or bake
  the tube. Oxygen and carbon radiation would dominate the energy loss at these
  temperatures.
- **Gas use:** under 1 mL STP per shot (static fill).
- **D2 at a few hundred eV produces some D-D neutrons.** This is for your
  radiation consultant.

## 11. What is still open before committing hardware
1. **Coil shape** is the one decision that is hard to change. The uniform coil
   (A) is cheapest and works with H2. The mild graded coil (B) costs about 1.7 ×
   the bank energy but works with both gases. Option A's bank, rated 20–25 kV,
   drives B.
2. **Pre-ionization physics** (frc-ibl.4): the ionization fraction and trapped
   flux against PI settings are empirical. Plan to scan them.
3. **Running now:**
   - Hall physics and the toroidal field;
   - the central-rod toroidal-field experiment.

   They will update §7 and §9.
4. **Lifetime** is not predicted by the 2D model. The n=2 mode and anomalous
   transport set it; the quadrupole is the main lever (§6). The 2D lever scan
   (EXPERIMENT_PLAYGROUND.md, Lifetime levers) found:
   - 0.3 T is shorter-lived (§2).
   - The slow mirror and the re-boost give no gain.
   - The 1 kV flux core holds the flux flat, but the FRC shortens to E ≈ 1.2.
   - Plain option A at 0.2 T stays the baseline: 1.11 mWb and E 2.6 at 40 µs.
