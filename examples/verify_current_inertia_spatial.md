# Nonlinear inertia axis audit

Checks exact engine/reference/protocol and field hashes, all five resolutions,
finite arrays and actual mesh coordinates. Reintegrates the **full** convective
electric curl in physical r at order 8, independently of the measurement's
reduced field and order-12 x quadrature. Replays the reported forcing error.

The inspected scientific figure shows both global convergence and persistent
first-annulus error, plus the actual finest first-annulus reference/computed
profiles. A successful archive audit does not pass the axis-consistency gate.
The record explicitly distinguishes B* forcing from the physical magnetic rate
after the finite-inertia elliptic solve. It requires investigation before moving
integration and cannot be advertised as a formation-accuracy result.
