# Pressure-only primal/dual energy evolution

Holds density and ions fixed while advancing g/U with the verified pressure
midpoint and a separate **conservative face-flux** update of primal total energy.
The exact physical ion pressure should remain 16 Pa in this split. Any local
change is the measured discretization defect, not assigned physical heating.
Global ion energy stays fixed to the paired energy/flux solve accuracy.

All nonlinear residuals, physical energies, positive species minima, conserved
proxy and logarithmic entropy are stored at every accepted step. Initial/final
states and volume conventions are lossless; time and mesh controls are distinct.
No bulk acceleration, Hall cross term, resistivity or formation is implied.
The uniform ion pressure and nonlinear iteration budget are explicit case inputs.
