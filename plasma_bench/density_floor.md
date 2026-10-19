# density_floor.py

## Intent
A **declared** gap closure for field-reversed theta-pinch formation in the
fixed-chamber solver.

**The problem it addresses.** Behind the magnetic piston, the rising external
field drives an E×B drift of about 6e4 m/s away from the glass wall.
- Nothing replenishes the wall cells, so in NUCTE-III their density halves
  every ~75 ns: 8e19 → 8e17 m⁻³ in 0.5 µs.
- The explicit time step shrinks as √ρ, so the run stalls.
- In ideal MHD that region becomes a vacuum gap, which a fixed fluid grid
  inside the tube cannot represent.

**What really fills the gap.** In the real tube it holds residual gas, plasma
ionized from neutrals at a rate of about n_n·S(T_e). That rate depends on the
gap's electron temperature, which this single-temperature model cannot
supply. The closure stands in for that ionization.

**What it is not.** This is not the forbidden quiet numerical atmosphere. The
exterior outside the tube stays exact vacuum. Inside the tube, the minimum
density, temperature and ion mass are explicit parameters; every added kg
and J is counted and excluded from the conservation residuals; and results
must be shown insensitive to the floor level.

## Equations (Lie split after each accepted MHD step)
Cells with ρ < ρ_min receive Δρ = ρ_min − ρ of ions and electrons at rest at
`floor_temperature_ev`.
- **Momentum** is unchanged.
- **Kinetic energy:** the mass-loading drop becomes heat in the same cell, so
  total energy is unchanged by it.
- **Thermal energy added:** 2 e T Δρ/(m_i (γ−1)).
- **Cold dual-energy cells** update ρA from their entropy thermal energy.
- **Reported:** `floor_mass_kg`, `floor_energy_j` and `floor_cell_steps`.
  `mass_relative_change`, `fluid_energy_balance_residual_j` and
  `response_energy_balance_residual_j` exclude these counted additions.

## Verification (`python -m examples.dual_energy_verification --floor`)
- **Inactive floor:** bit-identical to running without a floor.
- **Active floor** (every cell raised, ρ 1 → 1.05): mass change excluding the
  counted floor mass is 5e-17, and the fluid energy residual is 2e-16 J.
