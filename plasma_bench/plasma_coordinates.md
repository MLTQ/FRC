# Physical active coordinates

`PlasmaCoordinates` removes only the already constrained transverse axis and
magnetic flux slots. Fixed scales use initial volume-mean number density,
physical sound speed from mean separate pressures and equal-energy magnetic
units. n/C scale with n0, regular q/l with n0 c/R, mz with n0 c, psi with
sqrt(mu0 M n0)c R², and g with sqrt(mu0 M n0)c/R. Bias stays separate.

`pack`/`unpack` are linear coordinates, not modified physical variables or
positivity floors. `background` supplies an explicitly approximate uniform
state for a solver preconditioner, never an accepted/evolved state. The original
spatial equations still act on the actual density and both actual temperatures.
`field` retains the existing fixed-flux harmonic representation. Scales remain
fixed for an integrator instance and are part of its residual contract.
