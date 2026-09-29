# Resistive component verification

Independent spherical continuum fields, current and cell heat must converge
under refinement, separately from exact discrete magnetic-loss/heat balance.
Physical B² plus analytic exterior energy independently checks source work on
a nonuniformly deformed mesh. Zero eta preserves even tiny material coefficients
under a much larger bias change, without cancellation or heat.

A separately diagonalized fixed-geometry semidiscrete eigenmode supplies exact
time integration targets; this isolates second-order midpoint time error and
is not a continuum spatial-accuracy claim. Invalid inputs leave input arrays
unchanged. No magnetic monotonicity or arbitrary-timestep accuracy is inferred
from energy dissipation.
