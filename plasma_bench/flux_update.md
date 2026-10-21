# Shared conservative/CT update

`FluxUpdate` stores one flux per radial/axial face, one electric field per
corner and the cylindrical radial-momentum source. `rhs(mesh)` takes geometric
flux differences and returns conservative derivatives, corner-flux derivatives
and the outward total-energy rate. It is the sole assembly path for both full
orders and mixed local updates; energy accounting uses the actually selected
boundary fluxes.

`local_fallback(low, troubled)` substitutes baseline fluxes on every face and
corner incident to any troubled cell, and its metric source. Each troubled cell
then has exactly the low-order Euler update. Shared entities have a single
choice, so adjacent cells cannot disagree on their exchange; CT remains a curl.
The returned affected mask includes diagonal neighbors sharing corner fields.
Untouched entities retain the high-order arrays exactly. No updated cell state
or energy is overwritten. The routine does not decide admissibility or silently
expand the region; `flux_limiting` owns that policy and its checks.

An optional `BoundaryPoynting` supplies impermeable-boundary energy fluxes from
the selected corner electric field and a fixed variational magnetic trace.
`rhs` inserts these actual face fluxes before taking the same conservative
divergence; it does not change a state or apply a cell-energy residual. Local
fallback keeps the boundary object and reevaluates its flux after selecting
corners. Both orders must refer to the same initial magnetic boundary load.

An annular mesh adds the inner-wall Poynting flux from
`BoundaryPoynting.energy_fluxes` on radial face 0. Its outward energy rate
includes −2π r_in dz Σ F_r(0), the energy leaving toward the axis.

`toroidal_axis` (2.5D, axis meshes only) is E_z on the axis. It enters the
B_θ/r update of the axis cells as −E_z(0)/∫r dr: the r = 0 term of
r F_r = −E_z, which r F_r cannot carry. The ideal and Hall parts vanish on the
axis; resistive η J_z(0) does not. `local_fallback` replaces it where an axis
cell is troubled.

`toroidal_boundary` (2.5D with a threading current only) holds the B_θ
Poynting flux on the wall, bottom, top and inner faces. `rhs` adds it to the
substituted variational poloidal boundary flux, so it reaches the outward
energy. It then also returns a fourth entry, the B_θ power into the gas. The
default return stays three entries. `local_fallback` replaces these faces like
the others.
