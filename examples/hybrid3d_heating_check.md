# hybrid3d_heating_check.py

Runs the actual kinetic-ion, thermal-electron stepper against an independent
force-free decay solution, with no imposed heating history or fixed ions.
`python -m examples.hybrid3d_heating_check` saves ignored comparison evidence.

In a fully periodic uniform plasma, B=A(t)(cos kz,sin kz,0) has J parallel B.
Ions remain at rest, electron heat stays uniform, and the continuum solution is
A′=−η(Te)k²A/μ0, Te=Te0+(A0²−A²)/(3μ0ne). An independent DOP853 scalar ODE
provides the reference, using continuum k rather than the grid's wave number.

Both constant η and temperature-dependent Spitzer η are run at 16, 32 and
64 cells per wavelength. Reported quantities are B and Te errors, total
energy error relative to initial magnetic energy, spurious ion kinetic
energy, cutoff-removed heat and div B. This also checks that resistive
electric work is not counted again as ion heating.

The fixture uses a periodic device with no wall or gap; it cannot validate
absorbing-boundary heat losses or establish FRC formation/confinement.

`entropy_contact(cells)` is an independent transport comparison: uniform
pressure, sinusoidal density/temperature and constant bulk drift should
translate exactly for one period, leaving pressure and velocity unchanged.
It measures entropy-advection diffusion and energy error in full kinetic runs.
The velocity reference is the prescribed constant drift, independent of the
mutable particle array passed into the stepper.

`decay(..., backend='metal')` and `entropy_contact(..., backend='metal')`
run the same independent physical comparisons through the opt-in mixed-
precision Metal transport. Default remains CPU. Small periodic fixtures
measure physical error, not representative GPU performance; actual device
startup replays are timed separately with synchronization included.
`gpu_moments=True` also exercises experimental atomic particle deposition,
separately from the transport-only path, in these same full physical replays.

`energy_advection` isolates the electron transport PDE at prescribed uniform
density and velocity. A Gaussian translates exactly around a periodic domain;
outer CFL≈1.5 requires energy substeps. It measures positivity, energy closure
and spatial convergence at32/64/128cells on both backends. This is an operator
verification fixture, not a force-balanced formation experiment.
All three comparisons accept `fused_arrays=True` to exercise the opt-in
float64 loop-fused CPU stencils against these same independent solutions.
