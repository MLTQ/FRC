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
