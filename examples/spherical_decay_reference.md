# Independent spherical magnetic-diffusion eigenmode

`SphericalDecay` gives a regular analytic axisymmetric magnetic field in a fixed
conducting sphere and its exact infinite vacuum dipole. Inside,
a=3 A j1(ks)/(ks) exp(-D k²t), k=n pi/R; outside a=a(R)(R/s)³.
The interface has continuous a and radial derivative, with no current sheet.
The central limit and derivative use a series to avoid cancellation.

`fields` returns SI flux per radian, Br/Bz, volume current Jtheta=k² r a/mu0
inside and zero current outside, plus Faraday Etheta everywhere. Vacuum E is
not inferred from an exterior Ohm law. `heat_density` integrates the analytic
local eta J² over finite nonnegative time. `integrals` gives independently derived full-space
magnetic energy, power and heat. The derivation is MOVING_RESISTIVE_PLAN.md.

The geometry is held fixed for a component induction test. Its magnetic force
is not balanced by a prescribed physical plasma pressure, so this is not a
free-plasma equilibrium or an FRC formation experiment. No numerical operator
or evolved solver data is used in constructing the reference.
