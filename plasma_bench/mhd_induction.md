# Constrained induction and physical resistivity
`ideal_emf` implements a four-state local-Lax–Friedrichs UCT corner electric
field with full directional jump dissipation. It does not use a simple average
of four fluid-face EMFs. `resistive_emf` computes scalar eta Jtheta from a
centered axisymmetric flux Laplacian, with zero axis and radial-wall EMF.
`resistive_energy_flux` supplies face Poynting flux to the conservative total
energy equation; no Joule source is double counted. `ohmic_power` independently
estimates physical heating with corner dual-volume quadrature; energy-partition
agreement needs refinement and is not enforced through a corrective source.
Axial arrays include the duplicate periodic endpoint. All outputs use SI.
