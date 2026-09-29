# Current/thermal protocol

Every case specifies physical masses, continuum reference, numerical mesh/time,
quadrature, nonlinear budget, density normalization and model explicitly. The
original nonbarotropic pressure fixture is unchanged (n=1e20 m^-3, R=.035 m,
pe from 5 eV, pi=16 Pa, original Btheta, zero poloidal bias). Both conservative
total-cell energy and separately evolved species energies are tracked.

Five mesh levels through nr128, time halving, massless joint controls, quadrature,
normalization and tolerance controls retain failures rather than skip them.
The joint massless control also reaches nr128. Strong waves are checked both at
the target density and at the explicitly labeled inertia-stress density.
Exact finite-amplitude waves cover both helicities, reversed axial propagation,
physical target density, stronger amplitude and a labeled n=1e18 inertia
stress fixture. The latter does not replace the original argon target. No
electron effective-mass adjustment is exposed. Time and mesh studies are distinct.
