# Fixed-current entropy-space refinement

## Purpose
Refine only the diagnostic fluid representation from 8 to 16 to 32 radial
elements, with the original B32/64 current fixed in every coefficient. This
measures how much of its entropy source the fluid space can represent; it is
not a new coupled evolution or an FRC formation calculation.

## Protocol
The original base is checked bitwise constant in radius before replicating its
exact coefficients to each new fluid grid. No polynomial transfer, fitted base,
current projection or magnetic change occurs. Both directions and species use
quadrature 16/20, actual physical masses, the original diagonal fluid volume,
and independent physical tensor-mass solves for separate diagnostic rates.

The strong source is fixed from the original F8 reconstruction for each order.
Report how the new base entropy gradients differ from that original target due
to floating reconstruction, rather than silently replacing it. Every current
fingerprint must equal the preceding mass audit's original reconstruction.
The F8 production and exact-mass coefficient vectors must reproduce the prior
audit bitwise. Full raw vectors and numerical metrics require a fresh replay.

## Contracts
`run` measures original-K, volume-advection/diagonal-mass, and exact-physical-mass
errors against the same fixed strong source. Conservation and analytic mass
matrix comparisons remain explicit. All field cells subdivide each fluid
element, so no integration cell crosses a fluid discontinuity. Only the scalar
species operator is refined; no coupled mesh-convergence or growth conclusion.
No unit tests, production source changes, accepted-state filtering or repairs.
