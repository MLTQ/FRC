# coil_circuit.py

## Intent
Turn required coil fields into buildable bank numbers, and drive the formation solver with
circuit-consistent coil currents:
- the inductance of single-turn sheet coils and wound coils;
- capacitor modules with stray inductance, start switches, crowbars, and bias and
  pre-ionization banks;
- axial sections and firing skew;
- the plasma's load on the banks.

It is a lumped, linear circuit. Its idealizations are stated below.

## Coil model (`Strap`, `section_coil`, `filament_inductance`)
**Geometry.** Coils are thin coaxial bands (`Strap`: radius, axial centre, width).
- Each band is split into equal-width filaments.
- Straps sharing a `section` and `turn` sit between the same collector plates.

**Inductances.**
- **Mutual, between filaments:** exactly 2π ψ(loop).
- **Same-radius pairs:** the logarithmic near field is replaced by its strip
  average, adding μ0 a ln(d/g), with g the exact geometric mean distance of two
  collinear strips.
- **Self:** μ0 a (ln(8a/g) − 2), with g = w e^{−3/2}.

**Equal flux.** On µs time scales the copper is a perfect conductor. Every
filament of one (section, turn) therefore links the same flux, and current
crowds toward the coil ends.

**Reduction.** With G the filament→turn incidence, H the turn→section incidence
(series turns) and M the filament matrix:
- turn matrix L_t = (Gᵀ M⁻¹ G)⁻¹;
- section matrix L = Hᵀ L_t H;
- filament current per unit section current D = M⁻¹ G L_t H.

Neighbouring sections redistribute each other's filament currents; each
section's total stays fixed.

**Outputs.** `SectionCoil.field(r, z)` gives per-section (ψ, Br, Bz) per ampere.
`sources(waveforms)` builds the filament `Coil`s for `VacuumField`.

## Circuit (`Module`, `Crowbar`, `simulate`)
**Elements.**
- **Module:** C charged to V, in series with stray L and R. It is switched onto a
  section's terminals at `fire_time`.
- **`unidirectional`:** the start switch opens at the first current zero
  (ignitron-like). After a crowbar this isolates the ringing capacitor.
- **Crowbar:** L_c and R_c across the section terminals. It closes `delay` after
  the capacitor voltage first crosses zero (≈ current peak), or at a set time.
- **One-way crowbar (`unidirectional=True`, ignitron/diode):** needs L_c > 0.
  - Once armed it conducts only freewheeling current. It releases when its
    current reverses, and re-closes when the section voltage forward-biases it
    again.
  - Without this, a later bank on the same coil (a re-boost) would be shorted
    by the first bank's crowbar.
  - A plug coil pushing flux back through the main coil would otherwise open
    the main crowbar for good.
- **Switch events:** start, crowbar arming, gate opening, release and re-close
  form one state machine. Every event whose quantity reaches zero together is
  applied together, and the transitions are logged (`crowbar_log`).

**Branch equations** (s = the branch's section):
- module: L_m İ_m + R_m I_m = V_C − V_s, with C V̇_C = −I_m;
- crowbar: L_c İ_c + R_c I_c = V_s;
- section: V_s = Σ L_ss' İ_s' + R_coil I_s + dΦ_p/dt, with
  I_s = Σ I_m − Σ I_c.

The system is integrated with DOP853 between switch events. Simultaneous events
(symmetric modules) are all applied at once.

`energy_residual` closes the initial capacitor energy against stored, dissipated
and plasma-delivered energy.

**Ideal-source rules.**
- Parallel ideal sources are singular, so module stray L must be > 0.
- More than one crowbar per section needs crowbar L > 0.

## Coupling to the formation solver (`CircuitStepper`, `LinearWaveform`)
The solver's coil boundary is a flux conserver. On its own that means zero
external inductance: the plasma's flux exclusion is answered instantly by an
image current.

The stepper makes the coil's flux follow the banks instead. Each accepted solver
step:
1. The plasma's flux linkage with each section is computed,
   Φ_p,s = Σ_j (2π/μ0)(Ku)_j ψ_s(j).
2. It enters as an image load, so that V_s = d(L I_coil + Φ_p)/dt with
   I_coil = I_vac + I_img.
   - Only sections with a closed circuit (`connected()`: a conducting module or
     a closed crowbar) answer a linkage change:
     dI_img,c = −L_cc⁻¹ dΦ_p,c (the driver's `coupler`).
   - An open section carries no net current. Its image is frozen, so the image
     and I_vac stay continuous when a section fires or opens.
   - The earlier forms both failed. Loading open sections algebraically was
     unstable. An offset formula I_img = −L_cc⁻¹(Φ − Φ_connect) made the
     already-connected images jump when a new section connected: the plug coil
     then fired and immediately opened.
3. RK4 advances the banks, with switch events located inside the step.
4. The per-section vacuum-current waveforms used by the next step are
   refreshed. They stay continuous: each keeps its own value at the step start
   and catches up with the circuit over one further step, with
   slope = circuit slope + (circuit − held)/step. A source that jumped at a step
   boundary would change the field with no matching work.

Power through the load, L İ_vac·I_img = −İ_vac·Φ_p, equals the solver's source
work exactly. The driver option is `drive="coupled"`.

## Verification (`examples/coil_circuit_verification.py`)

**Inductance.**
- The uniform-current filament model converges at second order to Nagaoka's
  exact sheet-solenoid inductance (relative error 6e-5 at 128 filaments for
  l/a = 8).
- N narrow series turns tend to N² × Nagaoka: 0.966, 0.977, 0.985, 0.991 at
  N = 5, 10, 20, 40.
- **NUCTE-III coil:** equal-flux 68.7 nH against the published 70 nH (Kuroyanagi
  2014). Wheeler's formula gives 73.5 nH. The equal-flux inductance is 3–5% below
  the uniform-current value.

**Circuits against closed forms.**
- LC and underdamped RLC: 5e-10 and 3e-10.
- Ideal crowbar at V_C = 0:
  - the coil decays as e^{−R t/L} to 1e-14;
  - the module rings as I_pk cos(t/√(L_s C)) to 2e-9;
  - the capacitor reversal −I_pk √(L_s/C) matches.
- Unidirectional switch: opens at t_cb + (π/2)√(L_s C), and the held voltage is
  exact.
- Energy residuals ≤ 6e-10.

**Sections.**
- Mirror-image halves driven identically reproduce the common-collector coil to
  6e-11.
- Skewed coupled halves close energy to 6e-10.

**Stepper.**
- Against `simulate` at identical times: 3e-10.
- With a prescribed image load, energy closes to 9e-12.

**NUCTE-III banks.**
- The published 32 kV, 67.5 µF fast bank gives a 0.64 T lossless peak at the
  coil centre. NUCTE rates it at 0.6 T.
- The published 2000 µF, 2.4 µH, 8 mΩ bias bank peaks 103 µs after firing. The
  published bias rise is 90 µs.
- The main stray (27 nH) and pre-ionization stray (78 nH) follow from the
  published 4 µs and 124 kHz, so they are consistency values, not checks.

## Idealizations
- **Lumped and linear:** no switch arc voltage or resistance model beyond R, and
  no skin-effect resistance.
- **The plasma enters only** through its flux linkage, and the solver's flux
  conserver sits at its fixed outer boundary rather than on the coil surface.
- **Unidirectional switches** open exactly at current zero.
- **Crowbars** fire at the capacitor-voltage zero plus a fixed delay.
