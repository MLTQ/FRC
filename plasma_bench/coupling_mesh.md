# Fluid–magnetic mesh coupling

`CouplingMesh` aligns a finite-volume fluid cylinder with a larger nodal field
grid. The exterior has magnetic degrees of freedom only: no fluid density,
pressure or conductivity arrays. The existing symmetric stiffness K enforces
K u=0 on exterior nodes, with zero induced flux on the finite distant boundary.
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
not an evolving free plasma surface, open gas boundary or infinite magnetic box.
