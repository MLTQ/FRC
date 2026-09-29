# Sampling on common initial material coordinates

`locate` independently inverts the initial radial-sector Q1 map by resolving
each (r,z) point into its two angular ray directions, then finding its graded
radial shell. It returns cell indices, shape functions and reference derivatives.
Points outside the initial material polygon fail; only roundoff-sized boundary
differences are clipped. This is a sampling tolerance, not a plasma-state floor.

`nodal` interpolates nodal data, including velocity at reference vertices.
`fields` reconstructs actual Q1 flux gradients in the deformed geometry and
weak current, with P0 density/pressure/cumulative heat per volume. Interior
quadrature also returns the actual cylindrical-volume ratio dV/dV0 from the
two geometry Jacobians and r/r0, for conservative common-support integration.
These points avoid the collapsed center Jacobian. Coordinates/velocity can
also be sampled on the regular axis with `nodal`, which needs no inverse Jacobian.

Unlike equal normalized-sector coordinates, equal physical reference positions
identify the same initial material location across different shell grading.
The comparison support must lie inside both initial polygons; independent
polynomial/deformed-map checks validate inversion and field reconstruction.
