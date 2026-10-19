# dual_energy.py

## Intent
Opt-in dual-energy pressure for the fixed-chamber conservative MHD solver. In
cold or low-β cells, pressure recovered from total energy is a tiny difference
of large numbers: the relative error is about (magnetic energy error)/(thermal
energy). In the NUCTE-III formation the piston evacuates wall cells to
β ≈ 1e-6, and the scheme fails there.

## Equations
- ρA with A = p/ρ^γ obeys ∂(ρA)/∂t + ∇·(ρA v) = (γ−1) ρ^(1−γ) q_heat. Here
  q_heat is the physical Ohmic heating ηJ², averaged to cells from the nodes.
  The solver's numerical interface heating is excluded by design.
- Transport is first-order Rusanov: each face uses the larger adjacent
  |v_n| + c_f, and reflecting faces carry zero flux. The step is a single
  forward-Euler step from the start of each time step, which keeps ρA positive
  under the CFL bound.

**Why not advect with the solver's own mass flux?** Upwinding by the sign of
that flux injected downstream entropy against the flow at strong shocks,
because Rusanov numerical diffusion can reverse the face mass flux. That
variant failed the standing-shock fixture and was rejected.

## Selection (`cold_cells`, `select`)
- **Cold:** internal energy from total energy < `threshold` × total energy, and
  the cell is not a compressive shock cell. The shock detector is
  D = −∇·v · min(dr, dz)/c_f > `shock_threshold`, with positive internal
  energy. Negative internal energy is always cold.
- **Warm:** ρA resynchronizes from the total-energy internal energy, so shock
  heating enters the entropy.
- **Reported:** cold-cell counts, the largest entropy-thermal/total ratio in
  cold cells, and "disagreements" (cold cells whose entropy thermal energy
  exceeds 10 × threshold of the total).

## Measured choices (`examples/dual_energy_verification.py`, 16 cells)
- **Threshold 1e-3 is too low.** In the low-β coil pulse, total-energy
  pollution in the wall cells reaches 3.5–4% of the cell energy, so those
  cells read as warm and receive spurious entropy jumps of 2e3–1e4. The
  compatible-energy options (variational boundary Poynting flux, staggered
  energy) do not remove this pollution.
- **Without the shock detector,** the standing-shock pressure L1 error rises
  from 2.585 to 3.035. A detector at 0.05 flags weak low-β waves; 0.2 does
  neither.
- **Reset mode** injects about 0.13% of the total energy in the low-β fixture.
  Pressure mode conserves fluid energy to 1e-16 J and reports a −0.004 J
  cold-cell discrepancy.
- **Low-β adiabatic invariant:** error 3%. This is the mismatch between
  first-order entropy transport and second-order density.
