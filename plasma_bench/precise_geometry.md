# Geometry before conversion to physical coefficients

This module uses Python Decimal with50 significant decimal digits for physical
coordinates and their arithmetic. `coordinates` preserves exact binary float
inputs or saved decimal strings. The precision context is local to each
decorated operation; it does not change the application's global Decimal context.

`widths` subtracts endpoints before conversion to float64. `fraction` divides
coordinate differences before conversion. `quadrature` constructs physical
points in decimal geometry, keeping even very narrow intervals distinguishable;
its integration weights are physical float64 lengths. `interpolate` locates
points and forms relative coordinates before evaluating float64 fields.

This is a finite-precision coordinate representation, not spatial smoothing or
an arbitrary-precision fluid solver. Physical N/U/momentum and force/mass solves
remain float64 with the existing factored species scales. Positive widths that
underflow on conversion explicitly reject. Differences smaller than50-digit
geometry can still fail; no endpoint is moved apart, merged or assigned a floor.
Saved geometry must use decimal strings, not rounded float arrays or pickle.
