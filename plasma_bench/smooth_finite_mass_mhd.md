# Smooth finite-mass time integrator

Inherits the verified explicit unsplit RK4, state partition, source-work
quadrature, diagnostics and atomic rejection from FiniteMassMHD. Overrides
only the complete derivative with shared smooth bulk transport. The explicit
transport quadrature defaults to five Gauss points per direction.
Bulk electron entropy defaults to smooth conservative reconstruction; the
explicit limited comparison isolates its accuracy/dissipation effects.

This is a separately selectable smooth reference, without the baseline's bulk
Rusanov diffusion or limiter. It is not a new default for shocks or formation;
time/space refinement and positive-state checks remain necessary. No baseline
source or measurement is relabeled. The class does not enable physical transport
or vacuum coupling and has no independently proved general stability bound.
