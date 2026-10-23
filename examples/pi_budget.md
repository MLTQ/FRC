# pi_budget.py

Step 1 of bead frc-ibl.4: a 0D pre-ionization (PI) energy and ionization budget
for the 12-inch design, followed by an assessment of PI options. The module
docstring says, for each number, whether it is exact, computed or declared.
Output goes to `outputs/pi_budget/pi_budget.json`.

## Result (2026-10-22)
**Inventory.** 1.95·10¹⁹ atoms sit under the 1.2 m coil for option A (H2 at
3.4 mTorr) and 2.01·10¹⁹ for option B (D2 at 3.5 mTorr). The dossier's 60 J PI
bank stores **19 eV per atom**.

**Minimum cost.** The thermochemical minimum per ion pair (half the dissociation
energy, atomic ionization and 1.5 Te at 2 eV) is 18.8–18.9 eV.
- Full ionization would therefore need 59–61 J delivered at **100% coupling with
  no radiation or excitation loss.** That is not achievable.
- Bounds on the ionized fraction α:

| fraction of stored PI energy delivered (η) | E_c = minimum | 30 eV | 60 eV | 120 eV |
|---|---|---|---|---|
| 1 | 1.00 | 0.64 | 0.32 | 0.16 |
| 0.5 | 0.51 | 0.32 | 0.16 | 0.08 |
| 0.25 | 0.26 | 0.16 | 0.08 | 0.04 |

  (Option A shown; option B is within 3%.) Only the minimum column is a strict
  bound. A computed cost per ion pair E_c needs the H2/D2 reaction network
  (frc-transport.6). The private MCCC archives that network would use are no
  longer in `/private/tmp/frc-mccc-personal`, which the system cleared on Oct 22.

**Benchmark anchor.** NUCTE-III stored about 1.33 kJ of PI energy
(11.25 µF at 15.4 kV, ringing at 124 kHz, up to 135 kA) for 3.7·10¹⁹ atoms:
**227 eV per atom**, 12× ours. Its PI field, about 0.126 T, exceeded its bias
of at most 0.09 T. Both numbers are the published maxima.

**Flux freezing does not need full ionization.**
- Electron–ion plus electron–H2 collisions (σ_m = 1.5·10⁻¹⁹ m², declared) give
  a skin time of 10 main quarter periods (30.5 µs) at α ≥ 0.4% (0.8–0.9% with
  σ_m doubled) for Te = 2–5 eV.
- At Te = 1 eV even a fully ionized column falls short (28 against 30.5 µs).
- So the bias stays frozen as long as Te ≥ about 1.5–2 eV at main fire, not
  because ionization is complete.

**The real risk is a mostly neutral fill at main fire.**
- **Energy is not the problem.** Ionizing the rest takes only 2–12% of the
  3.2 kJ main bank (E_c 19–120 eV).
- **Time is.** Electron density in the H2 fill e-folds in 1/(n_H2 K_iz), with
  K_iz computed from the repository's NIST-BEB model:
  - Te = 10 eV: 1.45 µs (1.4 snowplow times);
  - Te = 20 eV: 0.48 µs (0.5);
  - Te = 50 eV: 0.23 µs.
- The implosion therefore ionizes the remainder in time only if its sheath
  electrons are at about 20 eV or more. The formation runs assume a fully
  ionized fill, so this is unmodelled.

## PI options (12-inch coil, 150 kHz, 40 nH stray, 9.8 µF; energy ∝ B_PI²)
| option | stored energy | eV/atom | B_PI | bank | constraint and cost |
|---|---|---|---|---|---|
| **dossier PI** | 60 J | 19 | 0.028 T | 3.5 kV, 31 kA | under the 0.05 T bias. Partial ionization only, likely ≲ 10–30%. |
| **NUCTE-III per atom, rung above the bias** | 0.71–0.73 kJ | 227 | 0.096–0.098 T | 12.0–12.2 kV, 107–108 kA | B_PI ≈ 2 × bias, so the field rings through zero during the PI, as NUCTE-III's PI did. The bias bank is unchanged. Needs capacitors rated for ringing reversal at 12 kV. |
| NUCTE-III per atom, kept under a −0.1 T bias | same | 227 | same | same PI bank | the bias bank goes from 5.3 to ~20 kJ (∝ B²), and the main bank must reverse a larger bias. The benchmark shows more trapped flux with bias (0.6 → 2.3 mWb from −0.04 to −0.09 T). Not circuit-checked here. |
| lower the fill | 60 J | 227 | 0.028 T | dossier | needs 0.29 mTorr, about 12× lower. This changes the whole design point, so it is not viable on its own. |
| z-discharge PI | not bias-limited | – | B_θ, not B_z | separate bank | needs electrodes at the tube ends (hardware is the user's call). Not modelled. |

## Recommendation (for the build dossier)
- **Build the PI bank for about 0.7 kJ stored:** the 9.8 µF bank charged up to
  about 12 kV, ringing at 150 kHz, about 107 kA peak.
- **Commission it upward from the dossier's 60 J (3.5 kV)** and treat voltage
  and delay as bench scan parameters.
- **Allow B_PI above the bias, as NUCTE-III did.** Fire the main bank after the
  PI rings down.
- **Keep the PI-to-main delay short,** so that Te ≥ about 1.5–2 eV at main fire.
- **Watch for a neutral residue** on the interferometer and camera. If the
  plasma density at main fire stays well below the fill, more PI energy is the
  lever.

Not established here: the deposited fraction η, the real E_c, whether ringing
through zero traps any PI flux of the wrong sign, and how a partially ionized
fill changes formation. A formation run with neutrals (frc-ibl.4 step 3) needs
neutral physics the 2D solver does not have.
