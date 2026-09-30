# Measured physical regular-fluid/current coupling

The sixteen-case instantaneous protocol and its full bitwise/physical audit are
in outputs/regular_coupling_verification. Frozen source:
8c9339a6eaf0c88b88bffe108d5ca2c535a8bfbf33fa4bd11b0c0aaebba6c784.
The full suite passes **402 tests in 183.03 s**. Equations and the explicitly
incomplete moving energy contract are in REGULAR_COUPLING_NUMERICS.md.

## Independent force and near-axis tests

The arbitrary smooth physical force has three independent components with
exponential/trigonometric radial structure and periodic axial variation.
Density varies in both directions. Exact regular momentum rates come directly
from the force components and physical mass, without using the work adjoint.

| Degree-5 radial elements | Radial force error | Azimuthal force error | Axial force error |
| --- | ---: | ---: | ---: |
| 2 | 2.18390e-5 | 1.74003e-5 | 1.30216e-5 |
| 4 | 8.37542e-7 | 4.65750e-7 | 5.03528e-8 |
| 8 | 2.84584e-8 | 1.38152e-8 | 2.73814e-10 |
| 16 | 9.10329e-10 | 4.26004e-10 | 8.67456e-12 |

Errors are relative L2 in the declared nodal physical volume. First-element
active radial errors decrease 7.28289e-9,7.14016e-12,6.85703e-15,2.85775e-15;
axial errors there decrease 1.46983e-5,6.23901e-8,1.13326e-12,5.39526e-15.
The last values are roundoff/quadrature limited, not unlimited asymptotic order.
Degree-3 radial errors decrease 8.61238e-3 to 1.63678e-5 and are also retained.

Maximum force-work mismatch is 5.75e-16, axial force mismatch 2.12e-16 and
angular torque mismatch 4.90e-16. These balances are necessary but distinct
from the analytic accuracy measurements. No force normalization was fitted.

## Actual Ampere current and original induction

An analytic regular cylindrical flux/toroidal field gives independent current,
electric, magnetic-rate and current-rate targets. Its declared manufactured
forcing uses the original finite-inertia equation. The discrete solve uses
actual fluid continuity for n_t, never the analytic target in its place.

Degree-5 runs with 4/8/16 fluid radial elements use 8/16/32 field radial cells.
Magnetic-rate errors decrease 4.62633e-4,1.58115e-5,5.05377e-7. Finest
current/current-rate errors are 4.73980e-7/4.73006e-7. The two-element fluid
case also uses eight field cells, explaining the visible coarse plateau.

Actual Lorentz-force relative errors at the finest matched resolution are
7.44700e-8,2.06750e-7,1.04214e-6 for radial/azimuthal/axial components. The
first-element errors are 8.76049e-8,2.23107e-7,1.38035e-6. Thus agreement in
arbitrary-force tests alone is not mistaken for exact magnetic-force accuracy.

At fixed eight fluid elements, increasing field cells 8→16→32 reduces magnetic
rate error 4.62627e-4→1.58115e-5→5.42625e-7 while the actual continuity error
stays 7.64366e-6. Increasing fluid axial points from 40 to 80 reduces continuity
error to 1.35843e-7. Separate fluid and field errors therefore remain visible.
Quadrature orders 8/12/16, density scales 1e16/1e18/1e20 m^-3 and zero axial
density contrast controls all complete. These are verification parameters,
not a statement that every such plasma obeys this fluid model physically.

Maximum original induction residual is 9.87e-12; the ideal/Hall coupled solve
reaches 1.09e-11, both below the unchanged requested 1e-10. Full conditioning
and component residuals remain in each raw record. The manufactured induction
energy identity agrees within 3.55e-16. The actual current-energy density
pullback agrees with an independent complex derivative within 4.36e-16.

## Work, missing transport and independent integration

The ideal force/electric work cancels within 9.18e-16 and Hall work vanishes
within 7.59e-17 in scaled backward error. With actual changing density, the
incomplete fluid + magnetic + current energy rate remains nonzero. At the
finest matched resolution it is −0.001498461177 W; the derived still-unclosed
density exchange is −0.001498461105 W. Maximum discrepancy scaled by all
work magnitudes is 3.13e-14. This is evidence of an explicit missing transport
contract, not a claim of full total-energy conservation or physical heating.

An independent reconstruction uses polynomial roots for radial fluid fields,
complex FFT coefficients for axial fields, and separately differentiated
compatible basis polynomials for electromagnetic fields. It reads the saved
coefficients without calling the production sampler, force or energy routines.
Order-16 versus order-20 physical integration agrees within 2.42e-15 for
current energy, density exchange, Lorentz work and manufactured energy rates.
The coupled current-energy rate differs by up to 2.52e-14; a nearly cancelling
coupled magnetic-energy integral in the zero-contrast control has relative
quadrature difference 1.47e-12. The maximum production/independent difference
is 3.68e-13 for that same small magnetic integral. Every value remains saved.

All sixteen original-platform array/metric replays pass bitwise, excluding
elapsed wall time. All new local jobs are terminal. The six-panel report was
visually inspected; its explicit mesh ticks distinguish fluid elements from
field cells. Earlier frozen sources remain unchanged.

NEXT derive the common physical current-advection and drift-stress work,
including the real continuity/density pullback and cylindrical terms. Then
complete species pressure/entropy and enthalpy/Poynting transport, independent
parent-equation rates and unforced nonlinear evolutions. This instantaneous
component has not replaced a coupled production solver. Exact moving exterior
vacuum, physical transport/ionization, real ends/coils, converged formation,
3-D stability and confinement are still unproved.
