# Experiment playground

The formation simulator is set up to try several experiments on the 12" device,
not just one. An experiment is a set of `Nucte` fields (geometry, gas, banks,
coils). Named experiments live in `examples/experiments.py`.

```bash
.venv/bin/python -m examples.nucte_formation outputs/my_run preset=frc12_uniform_h2 end_s=30e-6
```

Explicit `field=value` arguments override the preset. Tuple fields take JSON,
for example `aux_coils='[{"name":"plugs","z_center":0.635,"pair":true,"bank":{"voltage":15e3}}]'`.
Summaries come from `examples/design_scan_summary.py RUN ...` (table and
diagnostic signatures) and `examples/nucte_compare.py RUN`.

## Building blocks

**Main coil and bias.**
- Fields: `coil_radius`, `end_radius`, `taper_start` (graded coil),
  `centre_straps`/`end_straps`, `strap_width`.
- Banks: `main_capacitance`, `main_modules`, `main_stray`,
  `crowbar_inductance`, `crowbar_decay_s`, and the `bias_*` bank.
- The main and bias charge voltages are not inputs. They are solved so the
  vacuum centre field is `bias_t` at t = 0 and peaks at `peak_t`.

**`Bank`.** One capacitor bank: total `capacitance`, charge `voltage` (the
sign sets the field direction), `stray` inductance, `resistance`,
`fire_time` (relative to the main bank), `modules`, and a crowbar.
- The crowbar is one-way by default (ignitron/diode): it freewheels the coil
  current, releases if its current reverses, and re-closes when forward-biased.
- Start switches open at their first current zero (`unidirectional`).

**`extra_main_banks`.** Further banks on the main coil, for example a later
re-boost.

**`AuxCoil`.** An extra coil on its own section and bank: `turns` series
turns filling `length` at `radius`, centred at `z_center`. `pair=True` adds a
series mirror image at −z_center.
- Uses: staggered coils, end plugs or mirror coils, compression coils, or a
  solenoid inside the inner tube.
- All coils are mutually coupled through the equal-flux inductance matrix.
  A crowbarred main coil keeps its flux linkage, so it answers an aux coil's
  flux by changing its own current. The plugs below lower the main current by
  ~12% this way.

**`inner_radius`.** A coaxial insulating tube on the axis; the gas fills only
the annulus.
- The tube is snapped to a field row: r_in = i0 r_w/(nr + i0).
- It is reflecting for gas and transparent to field. It is not a conductor,
  so a metal tube is not modelled.
- Coils inside it need radius ≤ r_in − 2 dr.
- Closed flux is the ψ pocket around the O-point that touches no wall (inner
  tube, outer wall, ends).

**Drive.** Use `drive="coupled"` for any experiment with extra banks: the
plasma then loads every connected section through its flux linkage.

## Presets
| Preset | What it is | Banks |
|---|---|---|
| `nucte_case1`, `nucte_case6` | NUCTE-III benchmark cases (docs/NUCTE_BENCHMARK.md) | published NUCTE banks |
| `frc12_uniform_h2` | **Option A**: uniform 16.5 cm coil, H2 3.4 mTorr, 0.2 T | main 8 × 6 µF at 11.6 kV (3.2 kJ) |
| `frc12_graded_d2` | **Option B**: graded 20 → 16.5 cm coil, D2 3.5 mTorr, 0.2 T | 15.2 kV (5.6 kJ) |
| `frc12_graded_h2`, `frc12_uniform_d2` | the other two gas/coil pairings | — |
| `frc12_reboost` | option A plus a second main bank fired at 15 µs | 48 µF at 6 kV (0.9 kJ), 8 modules |
| `frc12_end_plugs` | option A plus one-turn plug coils at z = ±0.635 m (series pair), fired at 3 µs | 20 µF at 15 kV (2.25 kJ), 60 nH |
| `frc12_inner_tube` | option A with a 3" OD insulating inner tube (control) | — |
| `frc12_flux_core` | inner tube plus a 40-turn, r = 2.5 cm, 0.8 m solenoid inside it, fired at 4 µs in the main-field direction | 300 µF at 2 kV (0.6 kJ), quarter period ~65 µs, 14.6 kA, 1.75 mWb |
| `frc12_core_punch` | inner tube plus a 10-turn, 0.3 m coil, fast, fired at 8 µs | 1.6 µF at 19 kV (0.29 kJ), 2 µs rise, 24.7 kA, 1 T, 1.9 mWb |
| `frc12_flux_core_1kv` | the flux core on a 1 kV charge: holds the trapped flux flat (Lifetime levers) | 300 µF at 1 kV (150 J), 5.9 kA peak |
| `frc12_core_reboost` | 2 kV flux core plus a main-coil re-boost holding the external field near the control | core as `frc12_flux_core`; 8 × 6 µF at 3 kV (216 J) at 15 µs, no crowbar |
| `frc12_mirror_bias` | option A plus a slow mirror pair (z = ±0.635 m, r 0.17 m, 4 turns each) fired with the bias; bias and peak rebalanced | mirror 200 µF at 3 kV (0.9 kJ); bias −4.4 kV (19 kJ); main 12.2 kV (3.6 kJ) |

## Why a flux core in the main-field direction adds closed flux
With an inner tube, the separatrix level ψ_s is set at the tube side of the
FRC, while the O-point value ψ_O is frozen into the plasma: in ideal MHD
∂ψ/∂t = 0 at the O-point, where B = 0. A long solenoid inside the tube adds a
nearly uniform δψ outside its radius. That δψ lifts ψ_s but not ψ_O.

The closed flux 2π(ψ_s − ψ_O) therefore grows by 2πδ, against resistive
decay at the O-point:
d(ψ_s − ψ_O)/dt = δ̇ − r η |j_θ(O)|.

Sustainment needs a loop voltage 2πδ̇ at least equal to the resistive one.
Option A loses 2 mWb with a ~66 µs e-fold, i.e. about 30 V. The preset core
supplies ~1.75 mWb over ~60 µs. The opposite polarity removes closed flux.

The crowbarred main coil returns the core's flux across its own cross-section,
so the external field falls by about Φ_core/(π r_coil²): ~0.02 T here.

## Verification of the new capabilities
- **Axis case unchanged:** `inner_radius = 0` reproduces a NUCTE formation
  snapshot of the previous commit bit for bit.
- **Circuit:** extra banks and aux coils go through the verified circuit model
  (`examples/coil_circuit_verification.py`). Stepper against integrator agrees
  to 3e-8 (re-boost) and 2e-5 (plugs).
- **Annulus** (`examples/annulus_verification.py`):
  - **Static equilibrium:** holds to round-off. Max speed is 6e-10 m/s against
    v_A = 9e4 m/s.
  - **Resistive flux penetration** into the inner vacuum, against an
    independent 1D solution. It is second-order convergent: max error
    2.9e-3 → 7.1e-4 → 1.8e-4 of the induced flux at 1 µs for 12/24/48 cells.
    The inner-tube field is within 9e-5, and energy closes to 2e-11. Table in
    `plasma_bench/conservative_coupled.md`.

## Demo runs (2026-10-20)
All runs use option A: uniform coil, H2 at 3.4 mTorr, 0.2 T, 24×144 cells,
coupled banks. φ is the trapped (closed) flux, l_s the separatrix length, r_s
the midplane separatrix radius, and T = Te + Ti inside the separatrix. "Formed"
means closed, shorter than 80% of the tube, and held for 3 µs. Mutual coupling
between coils is included throughout.

| Experiment | φ at 10 / 20 / 30 / 40 µs (mWb) | At 20 µs: l_s, r_s, n, T | Flux e-fold |
|---|---|---|---|
| Option A baseline (`frc12_uniform_h2`) | 2.02 / 1.76 / 1.45 / — | 83 cm, 10.9 cm, 5.5e20, 147 eV | 66 µs |
| Re-boost bank at 15 µs (`frc12_reboost`) | 2.02 / 1.79 / 1.48 / — | 67 cm, 10.5 cm, 7.0e20, 169 eV | 72 µs |
| End plugs at 3 µs (`frc12_end_plugs`) | 2.02 / 1.75 / 1.45 / — | 76 cm, 10.7 cm, 5.6e20, 149 eV | 66 µs |
| 3" inner tube, control (`frc12_inner_tube`) | 1.24 / 1.01 / 0.90 / 0.71 | 49 cm, 10.6 cm, 6.3e20, 177 eV | 60 µs |
| Inner tube + flux core (`frc12_flux_core`) | 1.42 / 1.69 / 1.80 / 1.97 | 52 cm, 11.4 cm, 6.2e20, 175 eV | grows |
| Inner tube + core punch at 8 µs (`frc12_core_punch`) | 1.26 / 1.01 / — / — | 28 cm, off midplane | — |

**What each experiment shows**
- **Re-boost.** The second bank adds ~35 kA to the main coil at 15 µs. It
  compresses the FRC (at 30 µs: 49 cm × 9.7 cm, 7.0e20 m⁻³, 150 eV, against
  67 cm × 10.4 cm, 5.4e20 m⁻³, 129 eV). External compression cannot add closed
  flux: φ stays within 2% of the baseline. The hotter plasma decays a little
  more slowly.
- **End plugs.** The crowbarred main coil conserves its flux linkage, so firing
  the plugs drops its current from about 150 to 86 kA. The midplane field is
  unchanged, while the field at z = 0.6 m rises from 0.14 to 0.25 T, a mirror of
  ~1.3 at the coil ends.
  - The FRC is unaffected until it bounces back axially. Its re-expansion at
    12 µs then peaks at 1.01 m instead of 1.19 m.
  - Trapped flux is unchanged. To shape the field inside the main coil's
    footprint, aux coils must fire before the main crowbar, or the main coil's
    own shape must change (option B).
- **Inner tube.** A 3" tube costs ~40% of the trapped flux during formation:
  1.24 against 2.02 mWb at 10 µs. The closed flux is set by the flux level at
  the tube side. From ~20 µs, the inner separatrix lies on the tube surface:
  closed surfaces graze the glass, where real plasma would cool and pick up
  impurities.
- **Flux core.** The 0.6 kJ solenoid inside the tube (main-field polarity,
  fired at 4 µs, ~30 V loop voltage) more than replaces the resistive loss.
  - φ grows from 1.35 to 1.97 mWb by 40 µs, 2.8 × the control, as the reasoning
    above predicts.
  - The FRC also gets shorter and fatter (26 cm × 12.0 cm at 40 µs).
    Elongation near 1 is tilt-unstable in 3D, which this 2D model cannot show.
  - The main coil current also falls through mutual coupling, so the external
    field weakens.
  - A 1 kV core holds the flux flat, and a re-boost restores the external
    field, but neither keeps the FRC long (Lifetime levers below).
- **Core punch.** A 1 T, 2 µs pulse from a short (0.3 m) core coil at 8 µs.
  - It pushes the inner separatrix out from 4.5 to 6.8 cm by 14 µs.
  - Its fringing field then splits the closed flux off the midplane by 16 µs,
    and φ ends no higher than the control.
  - A fast, short inner coil disrupts the FRC rather than adding to it. The
    flux must come from a long, slow core.

**3D views and comparison chart.** All six runs were repeated with saved states
every 0.1 µs (`snapshot_s=1e-7`, `outputs/fine_*`) for smooth playback. They
agree with the 1 µs runs above to within 0.3% in φ and one cell in length.

```bash
.venv/bin/python -m examples.formation_display outputs/formation_display "outputs/fine_uniform_h2::Option A (baseline)" "outputs/fine_end_plugs::End plugs" "outputs/fine_reboost::Re-boost bank" "outputs/fine_inner_tube::Inner tube (control)" "outputs/fine_flux_core::Inner tube + flux core" "outputs/fine_core_punch::Inner tube + core punch"
```

This builds an overview page and one 3D viewer with a run selector. Serve `outputs/` locally
(the `formation-display` launch entry) and open
`http://localhost:8791/formation_display/index.html`.

## Lifetime levers (2026-10-20)
What the 2D model says about holding the FRC longer.
- **Setup:** every run uses option A at 24×144 with coupled banks and saved
  states every 0.1 µs, and runs to 40 µs, so the whole flux decay is seen.
- **Baseline:** `outputs/life_uniform_h2_p020` extends `fine_uniform_h2` and
  matches it bit for bit up to 30 µs.
- **Not modelled:** the n = 2 rotational mode and the tilt mode. The
  quadrupole (~7% of B_e) is the main experimental lever against n = 2, but it
  needs 3D and is out of scope here.
- **Elongation:** E = l_s/(2 r_s). Values below about 2 (marked †) are
  tilt-prone in 3D, which this model cannot show.

| Experiment (preset) | Banks beyond option A (main 8 × 6 µF at 11.6 kV, 3.2 kJ; bias 2000 µF at −2.3 kV, 5.3 kJ) | φ at 10 / 20 / 30 / 40 µs (mWb) | l_s × r_s at 20 / 40 µs (cm) | E at 20 / 30 / 40 µs | Flux e-fold, fit 8–30 / 8–40 µs |
|---|---|---|---|---|---|
| Option A, 0.2 T (`frc12_uniform_h2`) | — | 2.02 / 1.75 / 1.45 / 1.11 | 83 × 10.9 / 52 × 9.9 | 3.8 / 3.2 / 2.6 | 68 / 54 µs |
| Option A, 0.3 T (`peak_t=0.3`) | main at 16.2 kV (6.3 kJ) | 1.93 / 1.76 / 1.25 / 0.24 | 53 × 9.5 / 17 × 3.0 | 2.8 / 1.9† / remnant | 56 / 21 µs |
| Slow mirror before the main bank (`frc12_mirror_bias`) | mirror pair 200 µF at 3 kV (0.9 kJ, 11.7 kA); bias −4.4 kV (19 kJ); main 12.2 kV (3.6 kJ) | 2.07 / 1.86 / 1.48 / 1.04 | 76 × 10.8 / 40 × 9.3 | 3.5 / 2.8 / 2.1 | 63 / 47 µs |
| 3" inner tube, control (`frc12_inner_tube`) | — | 1.24 / 1.01 / 0.90 / 0.71 | 49 × 10.6 / 34 × 9.7 | 2.3 / 2.1 / 1.7† | 68 / 61 µs |
| + flux core at 0.7 kV, switch held closed | core 300 µF at 0.7 kV (74 J) | 1.31 / 1.32 / 1.26 / 1.17 | 52 × 11.0 / 31 × 10.8 | 2.4 / 1.8† / 1.4† | 1.2 / 0.49 ms |
| + flux core at 1 kV (`frc12_flux_core_1kv`) | core 300 µF at 1 kV (150 J) | 1.34 / 1.41 / 1.39 / 1.33 | 52 × 11.1 / 27 × 11.1 | 2.3 / 1.8† / 1.2† | flat |
| + flux core at 2 kV (`frc12_flux_core`) | core 300 µF at 2 kV (0.6 kJ) | 1.42 / 1.69 / 1.80 / 1.98 | 52 × 11.4 / 24 × 11.9 | 2.3 / 1.5† / 1.0† | grows |
| + 2 kV core and re-boost (`frc12_core_reboost`) | core as above; re-boost 8 × 6 µF at 3 kV (216 J) at 15 µs | 1.42 / 1.68 / 1.79 / 1.97 | 49 × 11.3 / 23 × 11.9 | 2.2 / 1.4† / 1.0† | grows |

**What each lever shows**
- **0.3 T shortens the life.** Up to 20 µs it matches 0.2 T (1.76 mWb). The
  plasma is hotter (200–270 eV against 130–180 eV) and denser, so flux-limited
  conduction along open field lines removes 2.2× the energy: 4.5 against
  2.0 kJ by 30 µs.
  - The FRC then contracts, from 53 cm × 9.5 cm at 20 µs to 24 cm × 6.7 cm at
    34 µs, and its flux collapses.
  - The 130 µs e-fold in DESIGN_SCAN.md was fitted over 8–20 µs only.
  - The late collapse is under-resolved (r_s below ~11 cells), but the split
    from 0.2 T is already clear at 25–30 µs.
- **A slow mirror fired with the bias does not form a mirror.**
  - The 24 parallel main straps link equal flux. They carry counter-currents
    (up to −1.07 A per mirror ampere in the end straps) that spread the
    mirror's flux over the whole bore.
  - At the centre the pair gives 3.9e-6 T/A, against 5.1e-7 for bare loops. At
    11 kA it would cancel the −0.05 T bias to −0.007 T. The preset therefore
    raises the bias bank to −4.4 kV (19 kJ, against 5.3 kJ).
  - The main rise then pushes flux into the mirror coils and cuts their current
    from 11.7 to 4.9 kA. On axis at peak, B(0.6 m)/B(0) rises only from 0.72
    to 0.88, never above 1.
  - In the run the FRC forms shorter and hotter and bounces less: l_s swings
    over 58–101 cm, against 64–122 cm. φ is 6% higher at 20 µs, equal by
    30 µs, and lower at 40 µs.
  - The model's straps are perfect conductors. Real strap resistance (loop L/R
    of tens of µs) would partly relax the counter-currents during a 100 µs
    ramp.
- **A flux core of about 1 kV holds the trapped flux.**
  - The core supplies an initial loop voltage of about 1.23e-7 Wb/A × V/5.77 µH:
    15, 21 and 43 V at 0.7, 1 and 2 kV.
  - The control's resistive loss is about 1.24 mWb/61 µs ≈ 20 V. Thus 0.7 kV
    decays slowly, 1 kV holds φ at 1.33–1.41 mWb from 10 to 40 µs, and 2 kV
    grows it, as the reasoning above predicts.
  - **Switch:** at ≤ 1 kV the plasma's back-EMF on the core almost cancels its
    early current. At 0.7 kV with the default one-way start switch, the switch
    opened at 4.6 µs and the run matched the control (1.25 / 1.03 / 0.91 /
    0.72 mWb, `outputs/life_core_0p7kv`). Held closed, its current dipped to
    −0.12 kA and recovered. At 1 kV it came within 23 A of opening.
- **The re-boost restores the field but not the shape.**
  - It was sized by replaying the coupled circuit with the recorded plasma
    image currents.
  - It keeps the main current within +2.0 / −2.8 kA of the control over
    15–40 µs. The core alone is −6.1 kA by 40 µs.
  - The FRC is unchanged to within 1% in φ and one cell in length, so the
    weaker external field is not what shortens it.
  - A re-boost with its own crowbar mainly adds a crowbar in parallel, halving
    the crowbar resistance, so this one has none.
- **Why every inner-tube FRC gets short and fat.**
  - By 40 µs, open-line conduction has removed 2.6–2.9 kJ in every inner-tube
    run, with or without a core, against ~1.0 kJ of Ohmic heating.
  - Radial pressure balance fixes the pressure at ≈ B_e²/2μ0, so the plasma
    energy scales as B_e² r_s² l_s.
  - Holding the flux holds r_s (11 cm against 9.7 cm), so the length must carry
    the energy loss: 27 cm at 1 kV, against 34 cm for the control.
  - Flux sustainment alone does not sustain the FRC. That needs energy (and
    particle) confinement or input, which this model does not supply.

**Recommendation.**
- **Lifetime baseline:** plain option A at 0.2 T, with no inner tube:
  main 8 × 6 µF at 11.6 kV (3.2 kJ), bias 2000 µF at −2.3 kV (5.3 kJ).
  - It has the most flux up to 30 µs and still 1.11 mWb at 40 µs.
  - It keeps E ≥ 2.6 throughout.
  - Do not raise it to 0.3 T for lifetime.
- **Best flux sustainment simulated:** option A with the 3" inner tube and the
  1 kV flux core (`frc12_flux_core_1kv`).
  - **Core:** 40 turns at r 2.5 cm, 0.8 m long, inside the tube, in the
    main-field polarity. Self-inductance 4.8 µH, plus a ~1 µH lead.
  - **Core bank:** 300 µF at 1.0 kV (150 J), fired at 4 µs. Quarter period
    ~65 µs; 5.9 kA peak in the circuit model; crowbar at its voltage zero.
  - **Switch:** one that stays on through a brief current reversal at turn-on,
    for example solid-state with an antiparallel diode. A switch that opens at
    the first current zero is marginal at 1 kV.
  - **Tuning:** 0.7–1.0 kV brackets flat flux.
  - **Result:** φ held at 1.33–1.41 mWb to 40 µs, 20% above the open bore at
    40 µs. The cost is the inner tube (about 40% of the formation flux, and a
    separatrix on the glass from ~20 µs) and an FRC that shortens to E ≈ 1.2
    by 40 µs.
- **Not worth building for lifetime:** the slow mirror (no flux gain; 0.9 kJ
  plus 14 kJ of extra bias) and the main-coil re-boost (no change in shape).

```bash
.venv/bin/python -m examples.formation_display outputs/formation_display_lifetime "outputs/fine_inner_tube::Inner tube (control)" "outputs/life_core_0p7kv_hold::Core 0.7 kV" "outputs/life_core_1p0kv::Core 1.0 kV" "outputs/fine_flux_core::Core 2 kV" "outputs/life_core_reboost::Core 2 kV + re-boost"
.venv/bin/python -m examples.formation_display outputs/formation_display_lifetime_open "outputs/life_uniform_h2_p020::Option A, 0.2 T" "outputs/life_uniform_h2_p030::Option A, 0.3 T" "outputs/life_mirror::Slow mirror + raised bias"
```

The 0.7 kV held run used the `frc12_flux_core` preset with the core bank at
0.7 kV and `"unidirectional": false` in its `aux_coils` JSON. Open
`/formation_display_lifetime/index.html` and
`/formation_display_lifetime_open/index.html` on the `formation-display`
server.

## Caveats
- **2D axisymmetric:** all coils are coaxial loops, and an inner tube or
  solenoid is a perfect cylinder. The n = 2 rotational mode is not modelled;
  an inner tube may change it, but this model cannot say how.
- **End boundaries:** the induced flux is fixed two cells beyond the tube ends.
  The vacuum inside the inner tube connects to the source field there, a
  stronger end constraint than a real open tube.
- **The inner tube** is an ideal insulator. Wall loading, sputtering and the
  solenoid's own heating and forces are not modelled. Plasma contact with the
  tube is the main practical risk of that experiment.
- **Late-time decay** is set mostly by the declared flux-limited open-line
  conduction (`conduction_limit=0.1`) and by Spitzer + Chodura resistivity. It
  has no anomalous transport or particle end loss, so the e-folds compare the
  levers with one another, not with an experiment.
- **Coil straps** are perfect conductors with equal linked flux. That is
  accurate for the µs main rise, but only a limit for the ~100 µs bias-time
  coils.
