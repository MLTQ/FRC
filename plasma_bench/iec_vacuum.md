# iec_vacuum.py

## Purpose
First electrostatic component for the positive-grid IEC/POPS side experiment
(frc-cac). It computes a full 3D vacuum field around the open wire frame,
without particles, plasma, imposed focusing potential or frequency selection.
It cannot establish ion trapping, heating, fusion or optimal plasma parameters.

## Model and approximations
- Grounded spherical chamber of specified radius: an explicit reference
  geometry, not a claim about the user's unspecified chamber shape.
- All 60 vertices and 90 edges of a regular truncated icosahedron are retained.
  Wires have circular cross section and rounded, intersecting junctions.
  Cage diameter includes the wire radius; the feedthrough/support is omitted.
- Solve Laplace's equation outside metal with phi=V on the cage and phi=0
  on the chamber. The method of fundamental solutions places finite line and
  point sources **inside metal only**. These are a mathematical basis, not
  plasma charges or a numerical atmosphere. No source is softened in vacuum.
  Each junction has five interior point sources; dense cap-surface sampling
  addresses the measured few-percent cap error of the initial one-source fit.
- Direct line potentials/gradients are analytic. The grounded-sphere image
  term is integrated by Gauss quadrature. Its convergence is measured.
  A source at the sphere centre uses the finite image limit -1/R, without
  a singular image coordinate or softened charge.
- Unregularized weighted least squares fits the approximate wire boundary.
  Finite basis/surface sampling errors remain; independent off-fit boundary
  residuals, rank and refinement must accompany results. No symmetry averaging
  or sector constraint is applied to coefficients, potential or electric field.
  The specified regular geometry and uniform electrode voltage are symmetric.
- Vacuum potential scales linearly with voltage. This does **not** imply that
  a plasma's electron well, oscillations or losses scale linearly with bias.

## Components and contracts
- `truncated_cage`: constructs reference geometry from chamber/cage/wire sizes.
- `Cage.wire_distance`: distance to rounded wire axes, for metal exclusion.
- `Sources.kernel`: harmonic potential kernel and its negative gradient,
  in units per q/(4*pi*eps0). Evaluation on a source is rejected.
- `fit_conductor`, `solve_cage`: return `VacuumSolution` with fit residual/rank.
  Solve a nonzero voltage basis; linear scaling also covers an electrode at 0 V.
- `VacuumSolution.evaluate`: physical phi and optionally E=-grad(phi).
  Callers must keep evaluation points outside the metal and inside the chamber.
- `capacitance`: Gauss-law total cage charge divided by applied voltage;
  no artificial charge-neutrality subtraction is used.
- `cage_discretization(..., holdout=True)`: independent boundary locations
  for approximation assessment, never supplied to the coefficient fit.

`examples/iec_vacuum.py` checks against independent concentric-sphere Gauss-law
potential, field and capacitance, then reports actual reference cage results.
No kinetic plasma/electron closure, neutral collisions or time drive is here.

## Measured reference and limits
`outputs/iec_vacuum_reference_v4/`: 15 cm cage, 2 mm wire, 40 cm spherical wall.
At 2/4/6 panels per edge, common off-fit boundary RMS errors are
1.63/0.161/0.198%, maxima 17.0/1.04/1.39%. Local error is not monotonically
converged. Centre potential and capacitance change only 0.014/0.017% from
4 to 6 panels: phi(0)/V=0.866075, capacitance=11.090 pF. These are static
geometry quantities, not trapping depths or a validated electrode-loss model.
Image-quadrature doubling changes phi/V by <=4.5e-16.
The last source fit uses 2,460 coefficients and 19,740 boundary locations;
18,422 separate, common assessment points expose finite local errors.
Near-wire gradients and electron interception require further resolution before
quantitative kinetic losses are inferred. Centre E=0 alone is not a trap test;
the computed interior radial field has both signs. This static result does not
exclude driven orbit stabilization or an electron-created plasma well.
