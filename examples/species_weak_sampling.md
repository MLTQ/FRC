# Independent full-vector weak entropy reconstruction

## Purpose
Evaluate the integration-by-parts identity in SPECIES_INTERFACE_DIAGNOSTIC_DESIGN
for every fluid coefficient. Root-product radial polynomials and real Fourier
interpolation are independent of production samples, gradients and adjoints.
No dense 3072 by 3072 operator or evolution solve is formed.

## Components
- `fourier`: real periodic value/derivative matrices, retaining the Nyquist cosine.
- `radial_faces`: fluid traces from the owning field cell's fluid element;
  coincident left/right points remain distinct, including discontinuous tests.
- `WeakSampler`: tensor volume/gradient adjoints, independently integrated L2
  potential projection, and the entire weak entropy action and its decomposition.
- `faces`: one-sided outward covectors, split into radial interior/axis/wall
  and axial interior/periodic-seam partitions. Axis/wall are evaluated, not fixed.

## Contracts
The production quadrature projection and independent q16/20 projections are
supplied separately. Returned vectors are K*s before the physical -a/(e V)
species multiplier. Components are advective + half volume divergence - half
all-face flux + transfer residual; they reconstruct the independently integrated
skew action. The direct action remains a separate comparison, not a correction.
All physical weights use cylindrical volume or face area, in SI units.

The face sum is an integration-by-parts term, not a replacement numerical flux.
The projection transpose can have global support. No trace is averaged or zeroed.
The current arrays are independent physical reconstructions of the original weak
Ampere coefficients. Scalar entropy means S(C/n); quotient diagnostics remain
separate. No solver state or production operator is changed.
