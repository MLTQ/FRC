# Linear Hall reference coefficients

`HallParameters` specifies uniform singly charged carrier density, ion mass,
signed axial bias and prescribed scalar momentum resistivity, all in SI.
`two_fluid` retains physical electron mass, total species mass density,
reduced-mass current inertia and the exact center-of-mass Hall coefficient.
`hall` and `mhd` explicitly select massless-electron comparisons; MHD also omits
Hall. These are mathematical model comparisons, not automatically valid regimes.
Pair-mass input makes the Hall coefficient zero without discarding inertia.
No zero-density vacuum, coefficient floor, thermal closure or charge-state
inference is accepted. `alfven` is signed for propagation under bias reversal.

`hall_modes` consumes these immutable coefficients. `hall_cylinder` reconstructs
physical velocity using the same total mass. See LINEAR_HALL_NUMERICS for the
parent-species derivation, boundaries and limited linear scope.
