# Independent finite-mass compressible mode verification

All three frequency branches at physical, equal and zero electron masses are
substituted directly into the original linear ion/electron momentum, pressure
and continuity equations. Analytic Bessel/Fourier differentiation checks
Faraday and Ampere independently of the six-mode generator; fourth-order
finite differences check those derivatives on their component scales. The slow
mode's nearly curl-free E makes a tiny-curl finite-difference residual lose
digits, so that cancellation is not used as the physical oracle. The independently
eliminated dispersion and radial boundary traces are retained. These checks
test the continuum oracle rather than substituting a production RHS into itself.
Continuity/pressure residuals use sums of the original derivative-term norms;
pair-plasma pure shear modes have zero density/pressure perturbation, so dividing
only by those perturbations would test eigenvector roundoff against zero.
