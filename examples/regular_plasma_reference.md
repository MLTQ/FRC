# Independent coupled fixtures

`reference` selects the existing analytically differentiated finite-mass linear
Bessel/Fourier wave, or the exact nonlinear translating two-pressure contact.
Physical electron and argon ion masses are retained; these dimensions and
densities are numerical fixtures, not proposed laboratory operating points.

`exact` supplies n,v,pe,pi, full E, perturbation B and Ampere J plus independent
time derivatives. Contact B/J perturbations and acceleration vanish exactly.
The wave is a linear solution: nonlinear discrepancy must be controlled by
amplitude reduction, and errors use perturbation rather than background norms.

`fluid_state` directly samples the actual regular point variables and physical
species EOS, with the two inactive axis slots zero. `initialize` returns the
reference, regular fluid, state, coupled system and harmonic magnetic field.
Poloidal flux uses analytic corner values and annular Btheta/r uses independent
order-16 physical cell quadrature. There is no manufactured source or matching
discrete eigenmode fitted into evolution. Compatible field degree stays five.
