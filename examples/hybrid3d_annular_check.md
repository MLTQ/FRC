# hybrid3d_annular_check.py

Concrete independent comparisons for the new inner-tube 3D device geometry.
Run `python -m examples.hybrid3d_annular_check`; output is ignored simulation
evidence in `outputs/hybrid3d_annular_check/comparison.json`.

- `ring_projection`: a circular current filament lies in the declared gas
  annulus. Its exact elliptic-integral vector potential supplies solenoidal
  boundary-normal B. Gap B is compared against the analytic Biot–Savart field
  on three grids. It also measures the error when the annular circulation is
  discarded, gap divergence, and curl of the discrete harmonic.
- `ballistic_wall`: a complete cold-ion hybrid run with E=B=0. Inward ions
  hit the inner absorbing wall at the independently known time (r0−a)/v.
  The measured impact-time error resolves staircase boundary error; number
  and absorbed kinetic energy are compared with their exact totals.
  The ballistic step is dx/(4v), so time sampling refines with the boundary.
- `driven_vacuum`: a complete empty-chamber run under a tilted uniform
  magnetic-field ramp, compared with exact B(t) in tesla and div B. No plasma
  or numerical atmosphere is present; there are exactly zero particles.

These are verification fixtures. They do not establish FRC equilibrium,
stability, growth rates or confinement, and do not replace configuration runs.
