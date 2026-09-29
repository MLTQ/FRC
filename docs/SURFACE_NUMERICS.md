# Surface resolution and comparison integration

`MovingDomain(..., radial_grading=p)` retains the same Q1 magnetic/material
method and introduces one explicit initial-mesh control. For normalized radial
index xi, material radii are R[1-(1-xi)^p] and exterior radial offsets are
(Rmatch-R)xi^p. Thus p>1 resolves both sides of the material-vacuum interface
more finely while coarsening the interior. The default p=1 uses the original
uniform coordinates exactly. Geometry validation rejects nonfinite/nonpositive
p and shell spacing that collapses at floating-point precision.

Pressure, density, flux labels, magnetic vacuum, physical eta and heat equations
are unchanged. Angular cells and the exactly curved distant matching sphere
remain independent controls. This is initial mesh grading, not remapping or
automatic adaptivity. Smaller cells can increase time and splitting sensitivity;
the physical timestep and diffusion-error tolerance retain their distinct roles.

## Compare the same initial material locations

The preceding resistive report used coarse-vertex mass-weighted velocity,
common normalized material-sector coordinates for field comparisons, and
conservative fine-to-coarse cell restriction for thermal comparisons. Those are
still the documented historical measurements. For different radial grading we
now evaluate solutions at common physical initial coordinates X0=(r0,z0).

Within an initial angular sector, write X0=c0 e0+c1 e1 in its two unit ray
directions. The nominal shell radius is s=c0+c1 and angular Q1 coordinate is
tau=c1/s. Locating s between the recorded radial shells gives xi. Q1 shapes
then evaluate the evolved material map x(X0,t), velocity, total magnetic
coefficient and weak current. Actual B uses the gradient in evolved coordinates.
Pressure, density and cumulative heat density now remain their true P0 cell
values on each mesh. This stronger pointwise comparison includes subcell
differences that conservative averaging intentionally removes.
This inversion is independently checked with known affine maps and polynomial
field/current data on unrelated graded meshes.

Different meshes have different cell edges. Simply using one mesh's Gauss nodes
can miss discontinuous pressure/density or magnetic-gradient pieces. We instead
intersect every initial coarse/fine convex cell pair, reject only disjoint
bounding boxes, triangulate each positive-area overlap and use positive Gauss/
Duffy quadrature. Exact clipped-polygon moments, a discontinuous indicator and
independent revolved-frustum volume validate the common integration mesh.

Initial mass weights are dM=rho0 2pi r0 dA0, with the fixture's rho0=1. Evolved
test-volume weights multiply dV0 by

```
dV_test/dV0 = (r_test/r0) det(dx_test/d(xi,tau)) / det(dX0/d(xi,tau)).
```

The integral must recover the complete test polygon's physical initial volume.
A coarser angular polygon is inscribed inside a finer one; the omitted fraction
of the finer initial volume is reported. Separate surface-vertex velocity errors
retain boundary-motion uncertainty outside the scope of a bulk integral.
The fourth/sixth Gauss-order check quantifies remaining smooth integration error.

These are differences between numerical solutions on common material support,
not absolute errors against an exact dynamic solution. They are not directly
interchangeable with the historical vertex-weighted percentages. Differences
in density/pressure changes can be large relative to a small perturbation even
when absolute densities/pressures remain positive.

## Physics and provenance limits

All physical and numerical controls remain recorded. Independent gas expansion
and fixed-geometry decay assess component accuracy after grading; actual source
and resistive-relaxation cases expose the unresolved thermal surface. Positive
heat and a small energy residual do not prove an accurate current-carrying edge.
The [conditional surface-compatibility audit](SURFACE_COMPATIBILITY.md) explains
why the initialized ideal FRC plus constant-eta heating is not a known smooth
resistive free-boundary reference.

Production and shared measurement sources are frozen while cases run. Additional
typed inputs can use the same frozen measurement entry point; the extra time
control records its full configuration and a separately hashed reproduction
helper. Overlap integrals have complete input/source identities and can be reused
without changing their meaning when report formatting changes. Report/plot
hashes remain separate. No older dataset is relabeled with a current fingerprint.
