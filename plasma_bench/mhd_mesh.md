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
