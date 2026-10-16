# Coupled moving reaction moments

`ReactiveMaterial` composes source-half / moving-full / source-half. Each source
step gathers the current distinct material supports, evolves the actual existing
atomic moment source, and applies material_reaction_projection. The next moving
step therefore uses genuinely changed cell masses, velocities, thermal energies
and support blocks; none is reset to its previous value. Transport uses the
common Poisson field and compatible finite-mass material equations.

`MaterialReservoirs` retains cumulative produced photon and prescribed heat
amounts per transverse area. Each newly gathered source starts these incremental
reservoirs at zero; the ledger adds its actual production. Photons are not
silently advected with matter or claimed escaped. Binding energy is computed
from the current physical atomic populations. `totals` adds all reservoirs to
the independently measured matter and interior/exterior field energy.

The source's fixed velocity coordinates must match the material set. Source
projection returns that same explicit frame. An initially empty product receives
its source-derived moments and support, never an assigned velocity or seed.
`last_overlap` retains the most recent evolved source state for diagnosis if
projection fails. Negative thermal energy or a reaction-modified time restriction
rejects; no clipping, subcycling, hidden retry or shortened requested duration.

This is a first coupled operator, not production-ready formation physics.
Projection can be inadmissible, repeated overlapping meshes can proliferate
tiny cells, and local numerical velocity mixing must converge. Disconnected
pieces cannot interpenetrate as if they were independent fluids of one species.
Real-rate/charged moving trajectories, collisions, radial geometry, actual walls
and all original FRC validity gates remain required after independent limiting
cases. Full goal acceptance is unchanged.

The ledger separately accumulates observed field-energy changes caused by source
projection/coarsening. This diagnostic is NOT subtracted from full_energy_j_m2
and never repairs physical heat. `last_projection` exposes the actual last cut
counts and field changes; the two source halves both contribute to the ledger.

Two additional ledger fields record numerical projection thermal gain and the
L1 amount of local total-energy redistribution. They are diagnostics only,
never additional physical reservoirs or terms subtracted from full energy.
The precise projected-frame research alternative supplies them. The historical
transfer carries these fields unchanged and does not claim to measure them.
