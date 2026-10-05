# Separately derived full linear spectrum

Assembles every active fluid and magnetic degree of freedom from
StationaryLinearPlasma, not from individual nonlinear finite differences.
All branches are kept. Independently differentiates the actual physical energy
into fluid EOS/kinetic and magnetic plus Ampere^T I Ampere Hessians; Cholesky
coordinates expose energy symmetry and the full stationary spectrum.

The actual nonlinear driver must converge to this separate matrix under centered
mixed perturbations at epsilon controls. Full electric pressure response is
checked too. The failed full coordinate-difference attempt is archived separately.
Captured input later identified a localized near-axis toroidal magnetic
perturbation: the physical and coefficient residual norms disagree. The initial
pressure-only hypothesis (also in a frozen source comment) was disproved.
This separate spectrum does not fix that nonlinear solver rejection.

The imaginary-axis RK4 limit applies only to this bounded uniform stationary
mesh and its spectrum. It cannot establish a nonlinear moving-state CFL or a
long contact/formation trajectory. Full operator, positive metric, spectrum,
random direction and nonlinear comparison are returned for independent replay.
