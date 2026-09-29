# Finite-time entropy transport verification

Checks the sheared continuum entropy equation by complex-step differentiation
and its pressure-energy/proxy integrals against independent closed Bessel/radial
expressions. A direct root-product physical Galerkin matrix exponential supplies
an independent temporal target; three RK4 step sizes must show fourth-order
improvement. An extreme step must reject without modifying accepted entropy.

Two physical grids verify actual finite-time entropy-wave and near-axis spatial
refinement, with physical logarithmic entropy conservation reported separately.
These are prescribed-current/density transport checks, not self-consistent
magnetic or moving bulk trajectories.
