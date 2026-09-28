# Vacuum boundary verification

`compact_current` specifies a smooth positive toroidal current well inside the
unit chamber. `infinite_flux` independently integrates analytic filament-loop
flux over that current distribution using Gauss quadrature. The integration
weight is Jtheta dR dz, the current carried by each circular filament; no extra
2 pi or radial Jacobian multiplies it. Probes lie outside the current support,
so the integrand is nonsingular. Increasing quadrature order checks its error.

`finite_box_flux` solves K u=mu0 Jtheta dA on an aligned finite field domain with
zero induced outer flux. Probe values are read at actual grid nodes, without
interpolation. Compare separately against infinite-space flux while varying
outer extent and grid spacing: satisfying K u=0 in vacuum does not eliminate
outer-boundary truncation error. This check is independent of the time integrator.
With `exact_boundary=True`, independent loop quadrature supplies the nonzero
Dirichlet boundary instead. That verification-only option separates spatial
discretization error from finite-box truncation. The dynamic solver continues
to use the declared homogeneous induced boundary.
