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
- `snapshots.npz`: total flux, density and pressure every 1 µs.
