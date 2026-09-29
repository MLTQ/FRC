# Moving meridional geometry and magnetic quadrature

`MovingDomain` supplies material Q1 cells initially filling a polygonal spherical
cross-section. Its center is a shared node and its axis stays at r=0. Material
nodes move in both r and z; there are no exterior gas cells or masses. Exterior
magnetic mesh nodes follow a declared linear extension of interface displacement,
with a fixed outer matching sphere. Invalid Jacobians/radii raise without repair.
All material nodes must remain strictly inside that sphere. Regular bilinear
cell corners are checked as well as quadrature points; only the intentionally
collapsed center edge may have zero corner Jacobian. The curved matching ring
is checked at quadrature points, not by a global mesh-topology theorem.

The last magnetic annulus is curved exactly onto that sphere. Its field trial
space blends inner linear nodal data into exact angular multipoles at the outer
surface, preserving continuity over entire faces, not just matching vertices.
For a=psi/r², regular vacuum modes are `(R/s)^(l+2) P_l'(z/s)/sqrt(norm_l)`
with norm_l=2l(l+1)/(2l+1). The infinite weighted energy has diagonal
`(l+2) R³`; the physical magnetic energy beyond R has diagonal `l R³`.
They differ by the integration-by-parts boundary term, which cancels globally.

`elements` evaluates positive quadrature, basis gradients and geometry gradients;
`stiffness` assembles integral r³ grad(a)·grad(b) plus the infinite exterior.
`Elements.shape_gradient` is its analytic virtual-work derivative at fixed field
coefficients. `volumes` and `volume_gradients` provide exact-to-quadrature
cylindrical cell volumes and their pressure-force derivatives. Material pressure,
mass and thermodynamics belong to the separate moving-plasma solver.

The initial material surface is polygonal and converges to the analytic sphere
with angular refinement. Only the distant matching surface is exactly curved.
Quadrature, material/angular resolution, exterior spacing/modes and mesh
distortion errors all require independent checks. This geometry alone does not
claim a complete free-boundary MHD solver or a mesh-remapping capability.

Optional `radial_grading` clusters both material and vacuum shells at their
shared surface. Material radii use R[1-(1-xi)^p]; exterior offsets use
(Rmatch-R)xi^p. The finite positive exponent p defaults to 1, retaining the
original uniform shell coordinates exactly. p>1 resolves thin surface layers;
it changes the mesh, never physical resistivity, initial density or cell heat.
Nonfinite, nonpositive or floating-point-collapsed shells are rejected. Angular
resolution, small center/axis elements and outer matching accuracy remain
independent controls. Independent physical volume, force, expansion and magnetic
decay checks must accompany use of a graded mesh.
