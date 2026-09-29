# Physical coupled-wave initialization and norms

Integrates density and all conservative momenta using independent physical
Gauss quadrature. Initializes exact corner psi, annular g, nodal electron and
cell ion pressures with the explicit discrete energy convention. The reference
is linear; finite-amplitude nonlinear residuals require amplitude controls.

Norms compare actual full coupled state with the unforced continuum wave after
subtracting only its known uniform background. Every denominator is the wave
perturbation, never background pressure or magnetic energy. B/J and nodal pe
use independent physical interpolation; fluid density/ion pressure and radial/
axial velocity use the declared cell representation, with regular vtheta=r Omega.
These field errors include interpolation and initial representation error.
Independent quadrature-order controls remain required for archived studies.
