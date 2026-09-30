# Conservative momentum and the physical bulk kinetic metric

Stores no alternate density. Given the common immutable `CompatibleDensity`,
sample actual cell means (n vr, n vz, n r vtheta) with regular powers 1,0,2.
Physical velocity divides these reconstructed momenta by the same n; angular
momentum is divided by r at positive physical quadrature. Total particle mass
is mi+me. K=integral M n |v|²/2, not kinetic energy of cell-center velocity.

`energy_gradients` pulls back physical velocity and -M |v|²/2 through the actual
momentum/density reconstructions, including volume. `energy_rate` keeps the
changing-density term separate from momentum work. `momentum_adjoint` maps
physical vector work to the conservative coefficients. None is a force update.

`gradients` supplies interior physical velocity derivatives, cylindrical curl
and divergence, with regular angular geometry. It does not substitute the strong
divergence for the single conservative continuity law. Polynomial jumps and
their possible weak-curl treatment remain distinct from these interior values.
All maps use the same physical geometry; nonfinite states reject without repair.
