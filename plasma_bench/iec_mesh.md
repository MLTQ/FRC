# iec_mesh.py

## Purpose
Optional Gmsh 4.15.2 adapter for a full 3D tetrahedral gas domain. The 60 rounded
junctions and 90 circular wires are boolean-unioned, then removed from a grounded
sphere. No symmetry sectors or gas cells inside metal are used.

## Components and contracts
`build_mesh(path, ...)` writes ordinary NPZ arrays: nodes, linear tetrahedra,
boundary triangles labelled 1=cage, 2=wall, and geometry/resolution metadata.
`iec_field` consumes these arrays without Gmsh. `shell` and `sphere` geometries
exist only for independent Gauss-law and manufactured-solution verification.

## Limits
Spherical chamber and absent feedthrough are explicit reference assumptions.
Piecewise planar boundary triangles approximate curved steel; wire loss and
near-wire fields require mesh refinement. Gmsh surface/volume sizing is not a
claim that every tetrahedron has that diameter. Actual element sizes are saved.
