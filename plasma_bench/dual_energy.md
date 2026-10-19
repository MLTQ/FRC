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

## Shock-override guard (2026-10-19)
The detector keeps the total-energy pressure only where the total-energy
internal energy is positive and at least `shock_min_fraction` (0.01) of the
cell's magnetic energy.
- **Why magnetic energy:** kinetic energy is recovered exactly from mass and
  momentum, so cold hypersonic shocks (kinetic-dominated) stay on total energy.
  Magnetic energy carries the discretization error, so magnetically dominated
  compression stays cold.
- **Failure it fixes:** in NUCTE-III with pre-ionization only under the coil,
  axis cells at z = ±0.80 m held β ≈ 1e-4 floor plasma under magnetic
  compression. The plain detector put them on total energy, and they went
  negative at 1.14 µs.
- **Rejected alternative:** a guard of 0.01 of the *total* energy disabled the
  detector at the standing shock (pressure L1 3.03 instead of 2.59).
- **Effect of the magnetic guard:** the shock fixture is restored to 2.5866,
  and detector 0.05 no longer causes spurious low-β heating.

## High-β reservoir failure and the entropy criterion (2026-10-19)
**Failure.** With `criterion="total"`, a cell's cold flag is computed from
E − K − M, which is the quantity that drifts.
- A cell that went cold in a low-β phase stays cold after the field leaves,
  because its E − K − M is still negative.
- This is true even when the cell is hot and β ≫ 1. In NUCTE-III case 1 at
  2.2 µs, such cells sit at r = 7.5 cm, z = ±0.57 m: 4e19 m⁻³, about 12 keV
  from the entropy pressure, B ≈ 0.02 T.
- In pressure mode the face fluxes of such a cell come from its entropy
  pressure. Its conserved E then pays the enthalpy outflow without limit
  (−6e8 W per cell), and the warm neighbours receive that energy.
- **Result at 3 µs** (16×128, Spitzer + Chodura):
  - cold-cell E − K − M is −14.8 kJ and warm cells hold +15.0 kJ;
  - the coils had delivered only 3.1 kJ;
  - the global budget still closes, so the response residual hid this.
- **Diagnosis by elimination** (single-step decomposition at 2.2 µs): the drift
  is unchanged with physical η off, interface η off, or with the LLF dissipation
  removed altogether. It is the central enthalpy flux of a cell whose E is
  decoupled from its pressure.

**Fix (opt-in): `criterion="entropy"` with `mode="reset"`.**
- **Cold test:** cold when the entropy thermal energy ρA ρ^(γ−1)/(γ−1) is below
  `threshold` × (that + K + M), and the cell is not a shock cell. A hot, weakly
  magnetized cell is therefore always warm.
- **Reset:** reset mode sets the conserved E of every cold cell to
  entropy thermal + K + M each step, so the switch back to warm is continuous.
- **Accounting:** the resets are counted as `dual_energy_change_j` and enter
  both the fluid and the response energy residuals as a declared source.

**Fixture check** (`examples/dual_energy_verification.py`, config `entropy_reset`):

| Fixture | entropy_reset | Pressure-mode default |
|---|---|---|
| Standing shock, ρ L1 | 0.39917 | 0.39915 |
| Standing shock, p L1 | 2.5867 | 2.5866 |
| Warm coil pulse | bit-identical to dual energy off | — |
| Low-β adiabatic-invariant error | 0.0305 | 0.0305 |
| Low-β reset correction | +0.004 J | — |
| Fluid residual | 5e-16 J | — |

**NUCTE-III case 1 at 3 µs:**

| Quantity | Fixed | Before |
|---|---|---|
| Physical thermal energy | 1.91 kJ | 19.9 kJ |
| Cold cells consistent with their entropy | yes | — |
| Reset corrections | net +193 J, \|Σ\| 681 J | — |
| Response residual, resets counted | 35 J (1.2% of the 2.8 kJ coil work) | — |

`select(..., criterion)` and `cold_cells(..., entropy, criterion)` take the new
arguments. The default stays `"total"`, which is bit-identical to before.
