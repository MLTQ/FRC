# Common physical radial momentum and kinetic metric

Interior radial faces store q=r n vr, with zero axis/wall faces constrained by
the space. The exact existing divergence-commuting particle-flux lift Qq is now
the physical radial momentum used for velocity, kinetic energy and force work.
Cell states retain actual means of n vz and n r vtheta. Radial cell means are
derived integrals, not an additional evolved or silently overwritten state.

`sampled_momentum`, `fields`, `energy_gradients` and `energy_rate` use the same
actual reconstructed density as thermal/current components and total mass
M=mi+me. The radial kinetic Hessian is positive definite on the independent
interior face DOFs. It avoids the rank-deficient lift of the old redundant
cell radial state. No momentum fitting or energy repair is applied to a step.

`gradients` differentiates Qq using dH=R D, including the −m_r/r term, and
differentiates the physical quotient by n. Other vector derivatives retain the
regular angular basis and periodic axial moments. `solve_rate` uses this true
variable-density kinetic metric to project physical radial particle-momentum
advection and subtract declared thermal/kinetic face work. Its physical Gram
preconditioner does not change density/mass; failed relative residual rejects.

This defines a new staggered radial state, not a reinterpretation of old saved
cell momenta. It still requires original-equation and full trajectory checks.
