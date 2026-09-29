# Curl-compatible poloidal energy

Interpolate psi linearly in x=r²/2 and z. Then Br=-psi_z/r and Bz=psi_x;
the interpolated field has zero pointwise divergence and matches all stored
face fluxes. Br uses shapes r_left/r (1-x_fraction), r_right/r x_fraction.
At the axis the left coefficient vanishes and Br is exactly proportional to r.
Bz is constant in x and linear in z. This is distinct from interpolating Br
linearly in r or squaring a cell-averaged field.

`energy_density` integrates this positive poloidal B²/(2 mu0) per physical
annulus. Radial mass coefficients use fixed 16-point Gauss quadrature of smooth
rational shapes (regular exactly in the first annulus); the linear-z squared
integral is exact. This fixed convention does not change with Hall quadrature.
`gradient` and `energy_rate` exactly differentiate that declared quadratic
energy. The duplicate periodic corner is folded once, not double counted.

The radial cache is bounded and read-only. This module does not alter a state
or repair its thermal energy. Any solver using this magnetic energy must
declare the choice consistently in construction, pressure recovery, total
energy flux and diagnostics. The old cell and staggered conventions persist.
