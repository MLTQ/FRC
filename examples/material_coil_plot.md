# Actual moving-plasma finite-coil results

Reads completed successful prescribed-pulse trajectories. Displays source
current, actual toroidal/meridional currents, heat/induced energy, both independent
work accumulations, and physical plasma Bz/Jtheta during the pulse. Spatial maps
are volume-averaged actual fields, with actual moving boundary and coil positions.
No fields are evolved, repaired or rescaled. Vacuum is not drawn as plasma.

At peak pulse, independently reconstruct induced flux at the physical coil
positions on the deformed vacuum mesh. Compare 2pi*sum(turns*psi) with the current
volume linkage integral used in source work. Their discrepancy measures a
physical finite-space reciprocity error; it never changes either simulation.
The report preserves actual summary/source identifiers and explicit limitations.
No formation, spatial convergence or physical transport claim follows from the
figure or integrated work agreement. Use same-grid time controls to assess the
inherited upwind time sensitivity before interpreting weak late-time currents.
For two otherwise identical inputs, actual saved states at quarter, half and
final duration give B/J/velocity material-label L2 differences, individual B/J
component errors (so dominant components cannot hide weaker ones), and maximum
position changes. No temporal interpolation or pass/fail tolerance is added.

For two different meshes with the same physical inputs, field/current/velocity
are evaluated in their actual basis at shared material labels. Ten-point Gauss
quadrature on each cell of the union of both reference meshes resolves their
piecewise maps. Volume weighting uses the final listed run's actual determinant.
Componentwise differences and maximum physical-position mismatch remain explicit;
two grids with unresolved time error do not establish spatial convergence.
The independent common quadrature is higher than the simulation's six points:
squared current norms include both the density profile and polynomial fields.

Integration method is included in each legend/summary. For split upstream runs,
the ambiguous instantaneous mutual/source-rate integral is absent, not zero;
its null value is omitted from the figure. The directly integrated external work
and every actual B/J/motion sample are still compared. Reconstruction needs only
the unchanged state/field definitions, not the historical time-step algorithm.
Different integrators on the same grid are labeled a method comparison rather
than a time comparison; agreement between them is not an analytic baseline.
