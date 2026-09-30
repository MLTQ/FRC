# Independent regular-coordinate checks

Exact monomial integrals verify Lobatto quadrature. Independent high-order
integration verifies that the positive-node first-element velocity basis has
exact physical x-weighted kinetic Gram matrix, with strictly positive eigenvalues;
there is no independent zero-energy axis velocity.

Original nonsingular material Euler equations verify all active state rates,
including strong near-axis convergence. Complex differentiation of independently
written physical H verifies local work; actual weighted angular momentum and
physical energy conservation are checked separately from accuracy.

The independently assembled radial linear equations must converge to Bessel
frequencies with a positive acoustic metric and skew physical work. Nonlinear
RHS difference controls, finite-time mixed runs and full continuous energy
audits are still separate requirements before adopting this prototype.
The nonlinear acoustic control perturbs both pressure/density and velocity;
matrix-exponential time refinement is checked. Nonzero inactive axis slots
reject explicitly, and a failed physical time step cannot mutate its input.
