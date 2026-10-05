# Matched original-continuum directional-rate verification

Samples the separately derived stationary contact's growing J0/J1 mode into
actual conservative fluid variables and analytic poloidal flux/cell toroidal
averages. Applies the declared wall-space *initial* approximation. Evolves no
fitted numerical eigenvector. Both positive and negative small perturbations
go through the actual complete moving finite-mass RHS; their centered physical
rate difference is compared with lambda times the original continuum fields.

Independent quotient/EOS reconstruction differentiates density, every velocity,
both specific entropies, B and current; full linear E is checked separately.
Reports each component's global/axis/wall/max rate and initial representation
error, quadrature16/20, mesh1/2/4/8/16 and amplitude1e-5/1e-4 control at mesh8.
No component is suppressed for a weak signal or large relative error. Raw states,
all physical rates/targets, integration weights and source/input hashes freeze.

The independent magnetic source has a measured small divergence residual;
poloidal flux initialization enforces divergence-free B with its own measured
representation error. Finite differences here test instantaneous consistency,
not complete perturbation trajectories or nonlinear stability. The bounded
contact fixture remains separate from formation and laboratory experiments.
