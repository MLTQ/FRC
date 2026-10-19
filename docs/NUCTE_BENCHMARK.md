# Benchmark experiment and solver map (frc-ibl.2, 2026-10-19)

## Chosen experiment: NUCTE-III (Nihon University)
NUCTE-III is a field-reversed theta pinch whose coil almost matches the 12"
design, and its formation measurements are published openly.

| | NUCTE-III | 12" design (FRC_DESIGN_WINDOW) |
|---|---|---|
| Coil radius | 0.17 m | 0.165 m |
| Coil length | 1.35 m between mirror points (18 + 2×5 straps of 5 cm) | 1.2 m |
| Tube | quartz, 0.23 m ID × 2.0 m | 0.30 m ID × ~2.2 m |
| Fill | D2 1.33 Pa (10 mTorr), base 1.5e-4 Pa | D2 ~3.5 mTorr |
| Bias | ≤ 0.09 T, ~90–110 µs rise | −0.04 to −0.06 T |
| Main | 0.6 T compression field, 4 µs rise; 67.5 µF / 32 kV bank, crowbarred | 0.2–0.3 T, 2–3 µs |
| Pre-ionization | ringing I_z discharge (2000); ringing θ PI at 124 kHz (2014) | ringing θ PI |

**Sources:**
- Takahashi, Iwata & Nogi, J. Plasma Fusion Res. SERIES 3 (2000) 498, "Source
  plasma formation of FRC for a long pulse operation".
- Kuroyanagi & Takahashi, Nihon Univ. CST proceedings O-10 (2014), "Analysis of
  circuit parameters for field-reversed theta pinch device".

### Measured targets
From Takahashi 2000, Table 1 and Fig. 2:

| Case | Mirror ratio | Field index | Axial contraction time (µs) | v_ax (cm/µs) | Reported flux Δφ (mWb) | r_s(0) (cm) | l_s (cm) | Lifetime (µs) | n=2 onset (µs) |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.10 | 0.03 | 6.5 | 3.8 | 0.15 | 6.2 | 62 | 65 | 28 |
| 6 | 1.80 | 0.30 | 3.0 | 8.5 | 0.4 | 8.0 | 38 | 50 | 22 |

Across all cases, n̄_e ≈ 2.2e21 m⁻³ and T_e+T_i ≈ 120 eV. For case 6:
- implosion with faster radial motion near the coil ends;
- reconnection complete by 3 µs;
- the contraction front reaches the midplane at 6 µs;
- equilibrium by 10 µs.

The diagnostics were an excluded-flux array, a side-on He-Ne interferometer and
visible light.

**Circuit (2014 paper):**
- Bias: 2000 µF, L = 2.4 µH, R = 8 mΩ, up to 60 kA.
- θ pre-ionization: 11.25 µF ringing at 124 kHz, up to 135 kA.
- Main: 67.5 µF with start and crowbar switches.
- Coil inductance 0.07 µH. Wheeler's formula gives 73.5 nH for this coil.

**Not published in these papers:**
- the bias field and main charge voltage used for each case;
- PI-to-main timing and crowbar timing;
- the definition of Δφ;
- the exact stepped-coil arrangement per case.

These are either recovered from fuller NUCTE publications or bracketed
by scans.

### Consistency with the design window
`frc_scaling` with nominal assumptions, the NUCTE geometry, D2 at 10 mTorr and
B_e = 0.34 T gives:
- **Density:** ⟨n⟩ = 2.29e21 m⁻³, against 2.2e21 measured. The measured shape
  implies particle retention f_N = 0.31–0.42; the nominal assumption is 0.4.
- **Temperature:** T_e+T_i = 118 eV. This agreement is partly circular, because
  B_e itself was inferred from NUCTE's n and T through the same ⟨β⟩ relation.
- **Classical τ_φ:** 72–115 µs, against measured lifetimes of 50–65 µs. This
  device class lives at roughly 0.5–0.9 of the classical estimate.
- **Trapped flux:** the rigid rotor (K=2) gives φ_p = 1.2–2.3 mWb, 3–10× the
  reported Δφ. Resolving this is part of the benchmark.

## Solver paths against the benchmark's needs
A field-reversed theta pinch starts with plasma filling the tube up to the
glass wall. The model needs:
- strong magnetosonic implosion shocks;
- reconnection at the coil ends;
- axial contraction;
- an insulating wall with vacuum between the wall and the coil.

| Path | Modules | Fit |
|---|---|---|
| **Fixed-chamber conservative MHD** | `conservative_coupled` (`ConservativeVacuumMHD`) | **Selected.** Eulerian r–z gas filling a reflecting, magnetically transparent tube. Shock-capturing finite volumes with positivity-preserving flux limiting and no floors. Analytic coils, and an elliptic magnetic-only exterior (fixed or multipole boundary). Single fluid and single temperature, with uniform resistivity, no Hall term and no toroidal field. |
| Split fixed-chamber MHD | `coupled` (`VacuumCoupledMHD`) | Same geometry with split hydro/magnetic steps. Superseded by the conservative path. |
| Finite cylinder (SBP, high-order) | `finite_cylinder_*`, `finite_coil_*` | Richer physics: two species, Hall, current inertia, Coulomb exchange, sealed ends. However, it has no shock capturing and has only run µT pulses. Not suitable for the implosion; a candidate later for post-formation equilibrium and decay. |
| Moving material (Lagrangian, free surface) | `material_*`, `material_experiment` | Full vector B/J, Hall, two temperatures. The plasma must be separated from the wall by vacuum, and runs stop at wall contact. A theta-pinch plasma starts at the wall, and the Lagrangian mesh cannot follow reconnection and contraction. Not suitable. |
| 1D slabs (Vlasov, reactive argon) | `vlasov_slab`, `reactive_material_*`, `species_material` | Component physics only; no coils or r–z geometry. |
| Stage 0–1, `AxisymmetricMHD` | `plasma`, `simulation`, `mhd` | Vacuum fields and diffusion, or verification-only MHD without vacuum coupling. |

## Gaps in the selected path, in the order the benchmark needs them
1. **Resolution and wall numerical diffusion.** The reversal current sheet sits
   at the wall, where the solver's declared interface diffusion acts. Measure
   with 16×128 → 32×256 → 64×512 runs.
2. **Resistivity.** It is uniform now. Add explicit Spitzer(T) from the single
   temperature, with any anomalous enhancement as a named, separately reported
   option.
3. **Ends.** Reflecting end walls retain particles that really stream out along
   open field lines. An end-loss model is needed before lifetime comparisons,
   but not before formation timing.
4. **Two-temperature physics and conduction.** These are needed for the
   equilibrium T_e/T_i split and the decay phase, not for the first formation
   comparison.
5. **Hall physics and toroidal field.** r_s/d_i ≈ 5 is marginal. It affects the
   reconnection rate and the toroidal field generated during formation.
6. **Pre-ionization (frc-ibl.4).** The initial state is assumed: uniform, fully
   ionized and frozen in the bias.
7. **Circuit (frc-ibl.5).** Coil currents are prescribed. The coil acts as a
   flux conserver through a fixed induced-flux boundary just inside it.
8. **Interferometer diagnostic.** Compute line-integrated density along the
   NUCTE chord for direct comparison.

## First coarse run (2026-10-19): reversal and early implosion, then a blocking failure
Driver: `examples/nucte_formation.py` with defaults: 16×128 cells, bias
−0.06 T, 0.4 T peak, uniform η = 1e-5 Ω·m.

**At 1 µs (241 steps, 13 s):**
- Field reversal at the wall extends over the full 1.41 m coil length.
- The excluded-flux and separatrix radius at the midplane is 10.9 cm (wall
  11.5 cm).
- The piston sheath has compressed to 1.5e21 m⁻³ while the wall side has
  emptied to 2.4e19 m⁻³.
- Mass is conserved exactly.
- **Numerical diffusion is large.** The response-energy residual is 8% of
  source work, and the solver's declared numerical interface heating is 29% of
  the coil's work (381 J of 1304 J). The reversal current sheet sits at the
  wall, where that diffusion acts. Results at this resolution are not
  quantitative.

**Failure at t = 1.614 µs** (3,111 further steps, dt collapsing to 75 ps):
- The first-order fallback proposal is itself inadmissible, and the solver
  correctly refuses to repair it.
- The failing cell is the outermost wall cell at z ≈ −0.43 m, where the
  mirror end straps begin.
- In that cell the piston has evacuated the plasma to 2.4e17 m⁻³ at 0.04 Pa
  under 0.79 T, so β ≈ 1e-6.
- Pressure is recovered as total minus kinetic minus magnetic energy, so it
  carries a relative error of roughly 1/β there. This is the known low-β
  positivity limit of total-energy MHD.
- A finer grid will not cure it: the evacuation is physical for a reflecting
  wall without a particle source.

**Conclusion of the map:**
- No existing path completes a field-reversed theta-pinch formation yet.
- The fixed-chamber conservative MHD path is the right one and runs the
  reversal and early implosion.
- Its blocking gap (frc-ibl.2.6) is a dual-energy (entropy) pressure update for
  low-β cells, verified on the existing shock fixtures.

**Possible remedies:**
- Dual-energy update. Its conservation consequence must be stated explicitly.
- Evolve gas energy, as the split solver does. That is immune to the
  subtraction, but it is the formulation the conservative path replaced.
- A physical particle source at the wall (ionization of neutrals and wall
  material). It is more realistic, but it is a separate model.

## Dual-energy fix and rerun (2026-10-19)

### What was added
`ConservativeVacuumMHD` gained an opt-in dual-energy pressure (`dual_energy=True`;
see `plasma_bench/dual_energy.md`):
- pressure mode, so total energy is never changed;
- cold threshold 0.1;
- shock detector 0.2.

It also gained opt-in timestep halving (`max_step_halvings`).

### Verification (`examples/dual_energy_verification.py`)
- **Off path:** bit-identical to the pre-change commit on the standing-shock
  and coil-pulse fixtures.
- **Standing shock:** pressure L1 2.5866 with dual energy, 2.5853 off.
- **Coil pulse:** identical to off.
- **Low-β coil pulse (β ≈ 2e-6):** off fails at t = 0.004. Dual energy
  completes, with fluid energy exact to 1e-16 J and an adiabatic-invariant
  error of 3%.

### NUCTE-III reruns
The low-β positivity failure is no longer what stops the run, but formation
still does not complete. Every variant fails with a first-order inadmissible
proposal:

| Variant | Fails at | Steps |
|---|---|---|
| Off: 16×128, η = 1e-5 (the original low-β wall cell) | 1.614 µs | ~3,350 |
| Dual energy, 16×128 | 1.277 µs | 755 |
| Dual energy + 6 step halvings | 1.277 µs | 755 |
| Dual energy, no interface diffusion | 0.786 µs | 52 |
| Dual energy, 32×256 | 0.686 µs | 187 |
| Dual energy, uniform η = 2.5e-4 (diagnostic) | 0.794 µs | 146 |

### Diagnosis of the dual-energy failure
- **Timestep.** The failing step is at about 1.24 µs with dt = 1.9e-11 s. The
  step is diffusion-limited by the declared interface numerical resistivity,
  which peaks at 0.082 Ω·m.
- **The accepted state is fine.** From it, both the high- and first-order
  single-stage updates of the failing cells change by only ~1e-4. Step halving
  to dt/64 does not help. The blow-up arises in the second stage, from the
  predictor state.
- **Current sheet.** Peak |J| is 2.75e10 A/m², a current sheet at the
  plasma–vacuum interface.
- **Failing cells.** They sit at r = 0.097–0.111 m, z = 0.63–0.68 m, under the
  mirror end straps. Their density is 2e19–6e19 m⁻³, and they are warm, not
  dual-energy cold. Their wave CFL is about 0.002.
- **What rules out simple causes.** Finer grids, removing the interface
  diffusion, and higher resistivity all fail earlier. That rules out simple
  under-resolution and the missing cold-edge resistivity as the cause.

### Leading hypothesis
The transparent-wall coupling is inconsistent under a fast external drive.
The CT electric field at the wall nodes, computed from face states, must match
the change the vacuum region imposes through the coil's flux. A mismatch
concentrates current at the interface nodes.

### Next action
frc-ibl.2.7:
1. Dump the predictor state of the failing step.
2. Locate the growing term.
3. Compare the wall-node E_θ with the vacuum −(1/r)∂ψ/∂t.

Then run the resolution study frc-ibl.2.1.

## frc-ibl.2.7 investigation (2026-10-19): root cause, fixes and wall E_θ check

### Root cause of the second-stage blow-up
The accepted state before step 753 was saved and replayed.
- The wall current is bounded (2–6e7 A/m²) for about 500 steps. Then, within
  two steps, the wall E_θ jumps from −5e4 to −7.9e8 V/m and then to 1.6e13 V/m.
- The second stage reconstructs the wall corner (16,11), at z = −0.83 m, with
  ρ = 1.14e-16 kg/m³, beside cells at about 1.1e-6. That is exactly the 1e-10
  relative trace margin of `mhd_positivity`.
- The resulting "fast speed" is 6.2e9 m/s, faster than light. The LLF UCT
  dissipation −½ a ΔB then gives E ≈ 7.8e7 V/m.
- The cell states themselves change by only about 0.1% from the accepted step.
  The CFL from cell averages cannot see trace speeds, and halving dt to dt/64
  does not help.

**Fix:** opt-in `local_trace_bounds`, which bounds reconstructed density and
thermal pressure below by the 3×3 neighbour minimum.
- Standing-shock accuracy is unchanged within 0.15%.
- On the coil pulse against a 64-cell reference, the density L1 error is 1.22%
  vs 0.74% at 16 cells and 0.49% vs 0.47% at 32: smooth-extremum clipping,
  vanishing with resolution.
- The step-753 blow-up disappears.

### Wall-node E_θ against the vacuum flux change (at step 752)
| Measure | Unbounded traces | Local bounds |
|---|---|---|
| Median \|E_wall − extrapolation from the vacuum nodes\| (fields ~1.55e4 V/m) | 88 V/m (0.6%) | 88 V/m |
| Midplane: E_wall vs first vacuum node vs coil alone | −1.60e4 vs −1.52e4 vs −9.1e3 V/m | same |
| Largest mismatch, just beyond the coil ends (\|z\| ≈ 0.78–0.83 m) | 9.9 kV/m, with sign flips | 3.9 kV/m, no sign flips |

- The difference from the coil-alone value is the plasma's flux compression.
- The wall coupling itself is consistent. It was corrupted locally by the
  unbounded corner traces.

### Second limit: gap evacuation stalls the time step
With bounded traces, the runs stall instead of blowing up.
- The rising external field drives an E×B drift of about 6.6e4 m/s away from
  the glass, in the wall cells beside the coil-end transition.
- With nothing to replenish them, their density halves every ~75 ns
  (8e19 → 8e17 m⁻³ in 0.5 µs), and dt ∝ √ρ.
- In ideal MHD this is a real vacuum gap between wall and piston, which a fixed
  fluid grid inside the tube cannot hold.

What physically fills the gap is ionization of residual neutrals, at a rate of
about n_n·S(T_e): roughly 2e6–8e6 s⁻¹ at 10–50 eV for half of a 10 mTorr fill.
That depends on the gap electron temperature, which this model lacks.
Charge-exchange friction cannot stop the drift: with low resistivity, the
frozen-in plasma moves at E×B, and the friction-to-slip factor
ηνρ/B² ≈ 5e-5.

**Closure:** `density_floor`, a declared minimum ionized density, added at rest
at the pre-ionization temperature. Its mass and energy are counted and excluded
from the residuals; the exterior stays exact vacuum. Runs use 1% and 5% of the
fill to measure sensitivity. A full neutral-ionization model is frc-ibl.4.

### Performance
`VacuumField` was 74% of the runtime: 16 evaluations per step, each looping
over 112 filaments. Grouping coils by waveform plus a time cache cut the step
cost from 75 to 14.7 ms, with fields equal to 1e-15.

## Formation campaign (2026-10-19, 16×128 unless noted) — SUPERSEDED
These runs used pressure-mode dual energy with the total-energy cold test.
The energy audit below shows that this combination created 4–5× the coil
work as spurious heat. The sequence and timing in the table may still be
qualitatively right, but their n, T, pressure balance, lifetime and energy
caveats are not valid. They are kept only for the record.
Common settings for every run in this section:
- dual energy, local trace bounds, step halving;
- D2 at 1.33 Pa, 2 eV pre-ionization;
- 0.4 T peak with a 4 µs rise;
- the declared density closure (1% of the fill unless noted).

Comparison quantities come from `examples/nucte_compare.py`.

| Configuration | End reconnection | Axial contraction | Outcome |
|---|---|---|---|
| Whole tube pre-ionized, uniform η = 1e-5, closure 5% / 1% | none for ~10 µs | none | Reversed flux fills the ~2 m tube. Density 2.0–2.7e21 (5%) vs 1.3–1.9e21 (1%); T 150–290 eV vs 420–580 eV: closure-sensitive. Flux e-folds in ~8 µs. |
| Whole tube, Spitzer η (± vacuum extended 31 cm past the tube ends) | none | none | The end boundary is not the cause. Temperature averages are polluted by hot, near-empty reversed tails. |
| Pre-ionized under the coil only (\|z\| < 0.7 m), uniform or Spitzer | weak / none | none | Even 1%-level plasma beyond the coil freezes the fringe bias flux. |
| Coil-only pre-ionization + Spitzer + **Chodura** (c = 0.1), case-1 coil (mirror 1.10) | ~3 µs, but at \|z\| ≈ 0.85 m, beyond the coil | none: the FRC lengthens to fill the tube | The FRC ends sit in the weak fringe field. |
| Same physics, **case-6 coil** (mirror 1.8; end straps at r = 12.7 cm, 12 filaments each) | ~3 µs at the mirror points | yes: 1.38 → 0.88 → 0.50 → 0.31 m at 3–6 µs | At 6 µs: r_s(0) = 9.0 cm, l_s = 31 cm. Then collapses by 10–11 µs; flux 2.4 → 0.9 mWb by 6 µs. Bias −0.06 T. |
| Case 6, **bias −0.09 T** (the published maximum) | closed by 2 µs | 1.25 → 0.44 m by 5 µs | 5–10 µs: r_s(0) 10.5 → 6.8 cm, l_s 44 → 22 cm, flux steady at 1.9–2.6 mWb, n 3.6–5.7e21. Collapses at 11–14 µs. |
| Case 6, peak field 0.5 T | as base | as base | Smaller, decays like the −0.06 T base. Bias, not peak field, controls size. |

**NUCTE-III case 6 (published):**
- reconnection complete by 3 µs;
- contraction front at the midplane at 6 µs;
- equilibrium by 10 µs, with r_s(0) = 8.0 cm and l_s = 38 cm;
- n̄ = 2.2e21 m⁻³, Te+Ti ≈ 120 eV;
- lifetime 50 µs, n=2 onset 22 µs.

**What matches at 16×128:**
- the sequence and timing: reconnection at the mirror points by 2–3 µs, axial contraction reaching the midplane at about 5–6 µs;
- the equilibrium size, for bias −0.09 T: r_s 7–10 cm and l_s 22–44 cm at 5–10 µs.

**What does not match:**
- density, 2–3× high;
- temperature, 3–4× high at 6 µs, which is the dynamic over-compression at contraction;
- lifetime, ~12 µs vs 50 µs.

### Numerical caveats that bound these results
1. **Wall numerical dissipation dominates the energy budget.**
   - The solver's interface numerical resistivity scales with the local signal
     speed, which reaches ~1e6 m/s in closure-density wall cells.
   - Its heat reaches 2× the coil's work by 12 µs in the bias −0.09 T run
     (12.7 kJ vs 5.5 kJ). It is 64% of the work by 3 µs with a 0.1% closure.
   - With it off, or with its speed capped at 2.4e5 m/s (the fill Alfvén
     speed), the formation phase is unchanged through 5–6 µs (r_s ≈ 10–11 cm,
     l_s ≈ 34 cm at 5 µs). The capped run's heat is ~29% of the work.
2. **End reconnection is partly numerical.**
   - At 32×256 (bias −0.06 and −0.09 T), the end reconnection has not
     happened by 3 µs, whereas at 16×128 it is complete by 2–3 µs.
   - The wall numerical heat falls to 5% of the work at 32×256.
   - The Chodura term (c = 0.1) does not yet dominate grid diffusion in the
     reconnection layer, so formation timing is not converged.
3. **The density closure matters.**
   - Particles added to keep the gap above 1% of the fill are 2–3% of the
     plasma mass.
   - Closure level changes density and temperature in the whole-tube runs.
   - In case 6 at bias −0.09 T, 0.1% and 1% agree at 2–3 µs.
4. **Physical resistivity inside the FRC is small:** about 1e-6 Ω·m (Spitzer
   plus Chodura, with v_d/c_s ≈ 0.05). The classical flux decay time is about
   0.4 ms. The observed decay is numerical or dynamic, not resistive.

## Energy audit, fixes and corrected formation runs (2026-10-19)

### 1. Spurious heating from the dual-energy cold test (frc-ibl.2.8)
A cell-by-cell partition of NUCTE-III case 1 (16×128, Spitzer + Chodura)
exposed the problem.

**At 3 µs:**
- Cold cells held E − K − M = −14.8 kJ and warm cells +15.0 kJ.
- The coils had delivered only 3.1 kJ, yet the plasma's physical thermal
  energy was 19.9 kJ.
- The global budget closed, so the response residual could not see this.

**Mechanism** (single-step decomposition at 2.2 µs):
- A hot, weakly magnetized cell had E − K − M < 0 and stayed "cold" by the
  total-energy test. One example: r = 7.5 cm, z = ±0.57 m, 4e19 m⁻³, about
  12 keV, B ≈ 0.02 T.
- In pressure mode its fluxes come from the entropy pressure. Its conserved E
  then paid the enthalpy outflow without limit (−6e8 W per cell), and the warm
  neighbours received that energy.
- Removing physical η, interface η or the LLF dissipation does not change this.

**Fix** (`dual_energy_criterion="entropy"` with reset mode; `plasma_bench/dual_energy.md`):
- The cold test uses the entropy thermal energy, so a hot cell is never cold.
- Cold-cell E is reset to the entropy value each step. The corrections are
  counted in both residuals.
- The standing-shock, coil-pulse and low-β fixtures are unchanged.
- Case 1 at 3 µs: physical thermal energy 1.91 kJ, down from 19.9 kJ.

### 2. Field/exterior energy closure: use the compatible path (frc-ibl.2.9)
On the historical cell / Riemann / RK2 path the response residual grows:
- **Case 6:** 1.4 × the coil work by 4 µs, roughly 10× smaller per grid
  doubling.
- **Case 1:** 15% of the coil work by 15 µs.

The repository's compatible path (staggered magnetic energy, variational wall
Poynting flux, RK3; `docs/COMPATIBLE_ENERGY.md`) closes it:
- 2.4e-4 for case 6 at 4–15 µs;
- 6e-5 for case 1 at 15 µs.

The graded-coil test (end filaments moved outside the field grid) changed
nothing on the historical path. That shows the error is the semi-discrete
boundary mismatch, not filament placement.

### 3. Remaining counted reset corrections (frc-ibl.2.10)
**Where they occur:**
- 88% of the reset corrections are in cells below 3× the floor density (the
  wall–piston gap). The FRC core gets none.
- They shrink about 3× per grid doubling.

**Effect of `magnetic_energy_flux="ct"`** (all magnetic energy transport as the
Poynting flux of the CT field, `vacuum_flux.md`): cold-cell drift falls 8×,
and 30× at the wall row.

**What remains:**
- Kinetic energy that discrete Maxwell-stress forces give the light gap gas
  without a matching E·J.
- First-order entropy diffusion from warm neighbours.
- Net size: +0.8 to +2.0 kJ by 15 µs in case 1 (floors 5% to 0.3%), against
  2.8 kJ of coil work.
- With open-line conduction on, the resets track the line loss (+3.6 vs
  5.0 kJ by 16 µs). The numerical gap heat is then carried out of the ends.

**The FRC is insensitive to all of this.** Case 1 at 15 µs across a 17×
range of floor density:

| Floor (fraction of fill) | r_s(0) | l_s | n̄ (m⁻³) | Te+Ti (eV) | φ (mWb) | Net reset (J) |
|---|---|---|---|---|---|---|
| 0.3% | 5.31 cm | 94 cm | 2.53e21 | 203 | 0.55 | +2035 |
| 1% | 5.10 cm | 100 cm | 2.54e21 | 198 | 0.52 | +1495 |
| 5% | 4.64 cm | 97 cm | 2.69e21 | 183 | 0.52 | +798 |

### Configuration used from here on
`examples/nucte_formation.py` with:
- Spitzer + Chodura (c = 0.1) resistivity;
- θ-PI under the coil only (|z| < 0.7 m);
- dual energy in reset mode with the entropy criterion;
- staggered quadrature, variational wall Poynting flux, RK3;
- `magnetic_energy_flux="ct"`;
- local trace bounds;
- 1% floor;
- 16×128 grid.

Comparison quantities come from `examples/nucte_compare.py`. Its flux e-folding
time is a least-squares fit of ln φ over the closed frames from 8 µs on.

### Results against NUCTE-III (Takahashi 2000)
**Case 6** (mirror 1.8, bias −0.09 T, assumed 0.4 T peak with a 4 µs rise):
- **Sequence:** reconnection at the mirror throats by 3 µs; the separatrix
  shortens 1.97 → 1.19 → 0.81 → 0.59 m at 2–5 µs, reaching about 0.5 m at
  6–7 µs.
- **Published sequence:** reconnection by 3 µs; contraction front at the
  midplane at 6 µs.
- **Then:** an axial bounce, with l_s back to 0.78 m at 9–10 µs and 0.50 m at
  15 µs.

| Quantity at 10 µs (equilibrium) | Simulation | NUCTE-III case 6 |
|---|---|---|
| r_s(0) | 7.2 cm (8.9 at 15 µs) | 8.0 cm |
| l_s | 78 cm (50 at 15 µs) | 38 cm |
| n̄ inside separatrix | 2.0e21 m⁻³ | 2.2e21 |
| Te+Ti | 320 eV | 120 eV |
| Trapped flux | 2.3 mWb (e-folding 47 µs over 8–15 µs) | Δφ 0.4 mWb (definition unclear) |
| Lifetime | (runs to 60 µs in progress) | 50 µs |

**Case 1** (mirror 1.10, bias −0.06 T): no axial contraction within 10 µs.
The separatrix fills the tube and shortens only from 11 µs, to 1.0 m by 15 µs.

| Quantity | Simulation | NUCTE-III case 1 |
|---|---|---|
| r_s(0) | 5.1–5.4 cm | 6.2 cm |
| n̄ | 2.0–2.5e21 m⁻³ | 2.2e21 |
| Te+Ti | ~200 eV | 120 eV |
| Flux e-folding | 38–42 µs | lifetime 65 µs |

**Peak-field bracket** (case 6, the drive amplitude per case is not published):

| Peak vacuum field | r_s(0) at 10 µs | l_s | n̄ (m⁻³) | Te+Ti (eV) |
|---|---|---|---|---|
| 0.40 T | 7.2 cm | 78 cm | 2.0e21 | 323 |
| 0.30 T | 7.6 cm | 84 cm | 1.6e21 | 253 |
| 0.25 T | 8.0 cm | 88 cm | 1.4e21 | 214 |

n̄T scales as B^1.7, i.e. pressure balance. No single drive amplitude matches
both n̄ and T: at NUCTE's n̄T (2.6e23 eV/m³) the simulated FRC holds about
40% fewer particles and runs about 1.8× hotter. The excess in T and l_s is
therefore a model discrepancy, not only the unknown waveform. Candidates, in
order:
1. Missing particle and energy loss physics: no radiation or ionization, a
   single temperature, no perpendicular conduction, and reflecting ends.
2. The assumed pre-ionization (fully ionized, 2 eV, under the coil only).
3. Coarse resolution: 32×256 is in progress.
4. The idealized mirror geometry: two coil radii instead of NUCTE's stepped
   elements.

**What is benchmarked and what is not:**
- **Robust:** the formation sequence and timing, r_s to 10–20%, n̄ to 10–20%,
  and a trapped-flux e-folding time of the same order as the measured
  lifetimes.
- **Not reproduced:** T (1.7–2.7× high) and the equilibrium length in case 6
  (2× long at 10 µs).
- **For design use:** take T from the NUCTE-anchored empirical scaling in
  `frc_scaling`, not from this model.
