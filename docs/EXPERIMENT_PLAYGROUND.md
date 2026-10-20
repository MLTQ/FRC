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
  - A smaller core bank (~1 kV), or a main-coil re-boost to hold the external
    field, should give flux sustainment without the shape change. This has not
    been run.
- **Core punch.** A 1 T, 2 µs pulse from a short (0.3 m) core coil at 8 µs.
  - It pushes the inner separatrix out from 4.5 to 6.8 cm by 14 µs.
  - Its fringing field then splits the closed flux off the midplane by 16 µs,
    and φ ends no higher than the control.
  - A fast, short inner coil disrupts the FRC rather than adding to it. The
    flux must come from a long, slow core.

The 3D views and the comparison chart come from `examples/formation_display.py`
(see its companion doc).

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
