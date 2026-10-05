# Admissible physical stage coordinates

Uses the frozen PlasmaCoordinates physical scales and removes outer-wall q.
An orthonormal null-space basis of the exact scaled toroidal trace parameterizes
every magnetic increment. Stages and endpoints are constructed in this space;
there is no post-step field/velocity repair. The original full coordinates are
retained for restricted preconditioners and diagnostic matrix comparisons.

`pack` is the coordinate pullback for tangent states/rates; `unpack` is the
admissible lift. Full field and fluid admissibility is separately validated by
the system and `field`, including floating trace tolerance. All nonlinear norms
are taken in this fixed physically scaled reduced coordinate basis. Tests must
also check full declared induction equations and discarded normal trace/work.
The inherited full scale/index metadata belongs to `full`, not reduced slots;
consumers must use `fluid_keep`, `magnetic_basis`, `nfluid` and `size` explicitly.
