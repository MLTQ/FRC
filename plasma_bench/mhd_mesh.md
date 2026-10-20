# Finite-volume MHD mesh
`MHDMesh` counts cells, unlike the Stage 0/1 nodal Grid. It supplies exact annular
volumes, radial measures and cell/face coordinates in SI. `face_fields` obtains
area-averaged Br/Bz from corner psi. `cell_fields` uses volume-aware face interpolation for primitive
states and energy; `divergence` computes the actual face divergence. The axis
has zero area and regular Br=0. Axial endpoints are duplicated periodic faces.
Cell Br uses the exact radial volume centroid in its linear face interpolation,
so a regular Br proportional to r is averaged correctly even in the axis cell.
No interpolation from the Stage 1 mesh is implicit; coupling needs an explicit
verified mapping. The MHD solver depends on these array shapes and measures.

## Annular chamber (`inner_radius > 0`)
Cells span `inner_radius..radius`, so `dr = (radius-inner_radius)/nr`. All
measures stay exact annular integrals: `radial_measure`, `radial_centroid` and
`volumes`. The first radial face is then a wall, not the axis. `face_fields`
divides by its nonzero radius, and `annular` reports the case. The default 0
reproduces the axis mesh exactly.
