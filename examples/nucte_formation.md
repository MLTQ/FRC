# nucte_formation.py

Coarse NUCTE-III field-reversed theta-pinch formation on the fixed-chamber
conservative MHD path (`ConservativeVacuumMHD`). This is the first frc-ibl.2
benchmark attempt, not a validated reproduction.

Run: `python -m examples.nucte_formation [OUTPUT] [field=value ...]`, e.g.
`... outputs/nucte_r32 nr=32 nz=256`. Default output: `outputs/nucte_formation`.

## Published inputs (Takahashi et al. 2000; Kuroyanagi & Takahashi 2014)
- **Tube:** quartz, 23 cm ID × 2.0 m, modeled as a reflecting, magnetically
  transparent chamber with reflecting end walls.
- **Fill:** D2 at 1.33 Pa, 300 K.
- **Coil:** eighteen single-turn 5 cm straps at r = 17 cm, plus five straps on
  each end at a smaller radius.
  - All straps sit in parallel on one collector, so the model gives them equal
    linked flux (solve `M I = 1`).
  - Each strap is split into four equal filaments.
  - The end-strap radius of 16.2 cm gives the case-1 mirror ratio of 1.10.
- **Main rise:** quarter period 4 µs.

## Assumptions (all explicit `Nucte` fields)
- Complete pre-ionization: uniform D+ at Te = Ti = 2 eV, current-free in the bias.
- Bias −0.06 T (the published maximum is 0.09 T).
- Vacuum centre field 0.4 T at the end of the rise. 0.6 T is the published
  compression field; NUCTE's pressure balance implies B_e ≈ 0.34 T.
- Crowbar L/R decay of 100 µs.
- Uniform resistivity 1e-5 Ω·m.
- Single fluid, single temperature, γ = 5/3, no Hall term, no toroidal field.
- Prescribed coil current with no circuit loading.
- Induced flux is held at zero on a fixed boundary just inside the coil
  (6 cells beyond the wall) and two cells beyond each end. This acts as a
  flux-conserving coil and conducting end planes.

## Solver options (Nucte fields)
- `dual_energy=True`: entropy pressure in cold cells, pressure mode, so energy
  is exact.
- `max_step_halvings=6`: timestep control.
- `interface_diffusion=True`: the solver's declared wall numerical
  resistivity. The run fails at 0.79 µs without it.
- `local_trace_bounds=True`: neighbour-minimum bounds on reconstructed
  density and pressure. Without them a wall-corner density of 1e-16 kg/m³
  produced a faster-than-light "fast speed" and a blow-up near 1.28 µs.
- `floor_fraction=0.01`: the declared gap closure, as a fraction of the fill
  density (see `plasma_bench/density_floor.md`).
- `resistivity_model`: "uniform" (`resistivity`) or "spitzer".
- `pi_half_length`: if > 0, the plasma is pre-ionized only for |z| below it,
  with a tanh edge of width `pi_edge_width`. This mimics a θ pre-ionizer
  under the coil; beyond it the plasma sits at the closure level.
- `end_cells`: vacuum cells beyond each tube end before the fixed
  induced-flux boundary.
- Boolean CLI values must be written `True` or `False`.

## Outputs
- `setup.json`: case parameters and the computed strap current distribution.
- `history.json`, written every 1 µs:
  - separatrix radius at the midplane, from the zero of the total flux off a
    reversed axis;
  - reversed and wall-separated lengths;
  - excluded-flux radius (r_w² − flux_loop/(π B_wall));
  - axis, maximum and minimum density, maximum pressure and speed;
  - mass change, response-energy residual, source work, Ohmic heat and the
    solver's declared numerical interface heat;
  - fluid energy residual, dual-energy cold cell-steps, cold-cell discrepancy,
    disagreements and step retries.
- `failure.json`: the accepted time and step count plus the error, if the run
  stops.
- `snapshots.npz` and `velocity.npz`: total flux, density, pressure and
  velocity, rewritten every 1 µs so runs can be inspected while in progress.

## Drives, extra coils and the inner tube (playground)
Named experiments live in `examples/experiments.py`. Run one with
`python -m examples.nucte_formation OUT preset=NAME [field=value ...]`;
explicit fields override the preset. See `docs/EXPERIMENT_PLAYGROUND.md`.

**`drive`**
- `"ideal"`: the prescribed centre-field waveform.
- `"circuit"`: capacitor banks through `plasma_bench/coil_circuit.py`.
- `"coupled"`: the same banks, loaded each step by the plasma's flux linkage.

Under `"circuit"`, the bias and main charge voltages are set on the main coil
alone, in vacuum: the centre field is `bias_t` at t = 0 and peaks at `peak_t`.

**Coupled loading.** Each section's image current changes as
dI_img = −L_cc⁻¹ dΦ_p, over the sections whose circuit is closed (a conducting
module or a closed crowbar).
- An open section's image is frozen. It therefore stays continuous when a
  section fires or opens.
- The vacuum current I_vac = branch − I_img stays continuous with it.

**`extra_main_banks`.** A tuple of `Bank`s on the main coil. Each has its own
capacitance, voltage, stray inductance, fire time, module count and crowbar
(one-way by default, so a later bank is not shorted by an earlier crowbar).

**`aux_coils`.** A tuple of `AuxCoil`s. Each is its own circuit section with
its own `Bank`: `turns` series turns filling `length` at `radius`, centred at
`z_center`. `pair=True` adds the mirror image at −z_center in series. A
negative bank voltage gives a field opposing the main field.

**Mutual coupling.** All coils are mutually coupled. A crowbarred main coil
keeps its flux linkage, so it answers an aux coil's flux by changing its own
current.

**`inner_radius`.** A coaxial insulating tube on the axis; the gas fills only
the annulus.
- The radius is snapped so the tube lies on a field row (`inner_rows`):
  r_in = i0 r_w/(nr + i0) with i0 ≥ 2. `setup.json` records the snapped radius.
- Coils inside the tube are `AuxCoil`s with radius ≤ r_in − 2 dr.

**Diagnostics.** An annular run reports `separatrix_mid_m` (outer separatrix),
`inner_separatrix_mid_m`, `trapped_flux_wb` and `inner_wall_bz_min_t`. Closed
flux is found by `nucte_compare.annular_closed_region`: the ψ pocket around the
O-point that touches no wall.
