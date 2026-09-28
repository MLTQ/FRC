# Topology diagnostics
`analyze_flux` finds spline critical points using gradient-bracketing cells,
classifies nondegenerate Hessians as O/X points, and detects physical axis Bz
zeros separately. It searches sampled contours enclosing O points for evidence
of closed flux. Counts are contour levels, not independent islands. Regions may
restrict detection to the prescribed plasma. Largest contours are lower-bound
geometric estimates, not automatically separatrices. Axis-attached psi=0 arcs
are closed to the axis only if resolved axis nulls support that interpretation;
otherwise separatrix estimates remain null. Off-axis X separatrix tracing and
robust treatment of degenerate/noisy/subgrid features remain research tasks.
Toroidal volume is 2pi times the polygon's radial first moment. Region crops use
the same roundoff-tolerant membership rule as the conductivity mask.
