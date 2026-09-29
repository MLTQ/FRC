# Independent Hall wave reference

Derives a continuum complex-frequency root directly from the species-reduced
quadratic, without importing solver matrices or coefficients. Stable root
selection avoids ion-branch cancellation. Both massless comparison models are
explicit. `CylindricalWave` uses independently evaluated J1 zeros, J1/J0 profiles,
and exact time dependence. Its curl is h K times the field; its velocity is
-k B0/(mu0 rho omega) times the field. Spatially varying positive eta J² heat
has an independent exact time integral, including the oscillating part.

This is an incompressible bounded-cylinder reference; the q=0 call to
`frequency` separately describes a Cartesian parallel transverse wave. Neither
is a formation experiment. `amplitude` has units T m (vector potential), unlike
the T amplitude of a magnetic perturbation. All modes are real fields obtained
from the complex phasor; helicity controls curl sign, not an ambiguous display
convention for right/left polarization. See LINEAR_HALL_NUMERICS.
