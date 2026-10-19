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
