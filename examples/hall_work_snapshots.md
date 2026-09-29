# Hall continuum cell snapshots

Tensor-Gauss integration independently averages the analytic HallManufactured
field and its exact instantaneous Faraday rate. g=Btheta/r uses annular volume
weighting; the r cancels in its integral. Density and B dot Bdot retain physical
r weighting. Corner flux and its rate are sampled directly. Defaults use order
10, separately controlled against order 14 in measurements.

This module imports no production curl, magnetic gradient, Hall projection or
Poynting formula. It defines references for genuine conservative cell averages,
not a reconstruction or an exact nonlinear time trajectory.

`wave_snapshot` similarly integrates a finite-pressure cylindrical mode;
`wave_hall_rate` differentiates only the linear Hall electric field by analytic
Bessel identities. It is distinct from the full mode time derivative, which
also includes ideal fluid induction. Comparison to the nonlinear operator
retains a declared O(amplitude²) residual. No reference matrix or numerical
current is used for the expected Hall Faraday rate.
