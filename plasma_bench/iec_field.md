# iec_field.py

## Purpose
Self-consistent electrostatics on a conforming full 3D tetrahedral vacuum/gas
domain. Linear finite elements solve -epsilon0 Laplacian(phi)=rho, cage phi=V(t),
wall phi=0. Real deposited charge is retained, including boundary shape weights.
No background charge, density floor, smoothing, symmetry filter or field clipping.

## Components and contracts
- `Field(mesh)`: assemble epsilon0 integral(grad Ni dot grad Nj), factor the free
  Dirichlet block once; construct exact face adjacency and conductor labels.
- `weights`, `deposit`: linear barycentric particle shapes, charge in coulombs.
- `solve(q,V)`: potential, constant-per-tetra E=-grad(phi), electrode charges
  (K phi-q) summed at cage/wall, field energy, space-charge energy.
- `locate`: initial/slice positions; nearest-centroid candidates with exhaustive
  tetra fallback. -1 means outside gas, not a silently relocated particle.
- `drift`: Numba straight-line face tracing to the first absorbing boundary;
  hit 1=cage, 2=wall, 3=unresolved traversal (caller must fail).

## Numerical limits
Linear FEM has piecewise constant E and finite particle shape/self-force error.
No ad hoc removal of a macroparticle's self-field is applied. Spatial, particle
and time convergence are required. Curved metal is approximated by triangles.
Sparse solves are float64 CPU; this is independent of the FRC hybrid/Metal path.
The field and particle shape use the same potential interpolation gradient.
