# Local complete-update admissibility

`admissible_update` checks density and internal energy after the **complete**
conservative/CT forward Euler proposal. A valid high-order update is returned
unchanged and does not evaluate first-order fluxes. Invalid proposals can raise
(`off`), use the historical global convex blend (`global`), or replace incident
face/corner fluxes and metric sources locally (`local`, the default).

For local fallback, independently assemble and check the entire low-order
proposal first. Flag every invalid high-order cell. `FluxUpdate.local_fallback`
chooses low-order values on all its incident entities, then recompute the full
mixed update. Add newly invalid neighbors and repeat. Previously flagged cells
use exactly the valid low-order update. Each unsuccessful pass adds at least one
cell, so at most the number of cells passes are needed; there is no iteration
cutoff or hidden global fallback. A first-order failure raises without repair.

This is binary local order reduction, not a claim of second order in troubled
regions or a universal MHD positivity theorem. It checks strict positive gas
internal energy and density, finite conserved values and cell magnetic field,
without absolute floors, pressure resets, entropy substitution or added heat.
The geometry/source/boundary baseline itself is conditional on the chosen time
step and physical state. All SSP Euler proposals and final state are validated.

Mass and total cell energy remain shared-flux conservative; corner CT preserves
compatible divergence. In a translationally invariant one-dimensional limit,
radial fluid and corner induction choices coincide. In genuinely 2-D mixed
stencils, positivity and conservation do not prove accuracy or preserve every
advected ratio; independent measurements remain required.

`LimiterReport` records troubled and incident-neighbor affected cell masks,
number of replacement passes (40 bisections for global comparison), and minimum
high-order weight. The solver accumulates these only after the entire step
validates. `invalid_cells` uses the same kinetic/magnetic energy convention as
`mhd_state.primitive`; it does not mutate state.
The explicit `energy_quadrature` selector is passed through every complete
proposal check. Both cell and staggered conventions are convex quadratic
magnetic energies of linear flux-derived fields, so the same SSP/convexity
argument applies. Selected boundary Poynting flux is reassembled from the
chosen corners by `FluxUpdate`, including all neighbor rechecks.
Magnetic units and quadrature reside in the shared energy helper.

`exempt` (optional cell mask) removes the internal-energy positivity requirement
for dual-energy cold cells, whose pressure comes from the advected entropy.
Density positivity and finiteness are still required everywhere. With
`exempt=None` the check is unchanged.
