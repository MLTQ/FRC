# Fluid–magnetic mesh coupling

`CouplingMesh` aligns a finite-volume fluid cylinder with a larger nodal field
grid. The exterior has magnetic degrees of freedom only: no fluid density,
pressure or conductivity arrays. The symmetric stiffness K enforces K u=0 on
exterior nodes. The default comparison fixes zero induced outer flux; the
multipole option below represents the field to infinity.
`harmonic_extension` preserves plasma-node values and solves that exterior.
At least two exterior intervals in each direction leave actual vacuum unknowns
between the chamber interface and the fixed distant boundary. Equal grid spacing
and aligned chamber end faces are required; unsupported geometries raise.

Node areas are intersections of node dual patches with the actual fluid region;
radial/end/corner interface nodes have half/quarter planar areas. `interpolate_r`
and `interpolate_z` map cell velocities to current-bearing nodes, impose zero
normal velocity at chamber walls and extrapolate tangential velocity linearly.
Their adjoints will distribute magnetic forces to cells with exact discrete
work pairing. `heat_to_cells` uses nonnegative patch fractions summing to one.

Field/current nodes and fluid cells are distinct. The interpolation's pointwise
axis/wall force error must be measured; energy pairing alone proves no force
accuracy. The interface is fixed and gas-filled, with reflecting ends. It is
not an evolving free plasma surface or open gas boundary.

`boundary="multipole"` replaces the finite outer condition with the variational
infinite exterior in `open_boundary`. The matching rectangle must enclose a
sphere containing the entire fluid chamber. `node_count` then includes modal
amplitudes after `grid_node_count` physical interior nodes. Modal rows have no
fluid mass, current density, heating, resistance or velocity interpolation.
Both boundary options satisfy the vacuum equation; only the multipole option
represents magnetic energy beyond the finite matching rectangle.
`interior/full` pack/expand the induced solution including modes; `grid_vector`
pads arbitrary physical grid data (notably external coil flux) with zero modes.
`harmonic_extension` preserves only plasma-node input values and solves vacuum
nodes and modes, ignoring supplied boundary data. The default fixed boundary is
retained as an explicit comparison. No arbitrary field may be packed as a
decaying induced trace; prescribed external bias/coils use `grid_vector`.
`complete_vacuum` works directly on unknown vectors, preserving active plasma
values and solving cached exterior/mode constraints. Its fixed factorization
and cross-block support fast midpoint iterations without changing the equations.

## Annular chamber (coaxial inner tube)
For `fluid.inner_radius > 0`, the field grid still starts on the axis. The inner
wall must lie on a field node row, `r_offset` = inner radius/dr, with at least two
vacuum intervals inside it; otherwise a ValueError is raised. The rows inside are
vacuum unknowns of K u = 0, like the exterior. The induced flux is fixed only on
the axis and the outer boundary.

The inner wall row is an active, current-bearing row with half dual area. It is
impermeable: no radial velocity is interpolated there. It is magnetically
transparent, with no conductor and no fixed flux. An insulating tube is
therefore modelled exactly. A metal tube, which would act as a flux conserver,
is not. `fluid_flux` returns rows `r_offset..r_offset+nr`.
