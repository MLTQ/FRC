# Constrained induction and physical resistivity
`ideal_emf` implements a four-state local-Lax–Friedrichs UCT corner electric
field with full directional jump dissipation. It does not use a simple average
of four fluid-face EMFs. Optional second-order reconstruction supplies both
quadrant fluid states and transverse staggered magnetic traces, retaining the
full one-dimensional Rusanov limit. `resistive_emf` computes scalar eta Jtheta from a
centered axisymmetric flux Laplacian, with zero axis and radial-wall EMF.
`resistive_energy_flux` supplies face Poynting flux to the conservative total
energy equation; no Joule source is double counted. `ohmic_power` independently
estimates physical heating with corner dual-volume quadrature; energy-partition
agreement needs refinement and is not enforced through a corrective source.
Axial arrays include the duplicate periodic endpoint. All outputs use SI.
`uct_electric` exposes the same four-state expression before boundary conditions;
the vacuum-coupled conservative solver applies its own reflecting-gas,
magnetically transparent boundary contract. `ideal_emf` retains its periodic
ends and perfectly conducting radial wall exactly.
`uct_electric(..., parts=True)` returns (central, dissipative) separately; their
sum is the unchanged UCT field.
