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

An explicit `FRC_GEOMETRY_BACKEND=mpfr` process environment selects the optional
`precision` dependency (gmpy2 2.3.2), with 192 binary bits, approximately 57
decimal digits. The default remains 50-digit standard-library Decimal. The
backend is fixed at import and cannot change inside a trajectory. Both local
contexts and scalar constructors request precision explicitly; MPFR's ordinary
53-bit default is unsuitable here. Invalid operations/division by zero raise.
Physical coefficients remain float64, and no spatial cut or population is removed.
`decimal` is the shared scalar constructor despite the historical name.

The complete matched 20 ns MPFR experiment retains identical saved physical
arrays and all 21,351 density cells, with coordinate differences around 1e-50 m.
All six full independent analytic cases now complete in `precise_material_mpfr`:
finest neutral/charged front errors 1.42e-15/6.94e-17 m and relative energy errors
5.68e-14/3.27e-15. These fixtures retain their disclosed 8-ULP cleanup option;
the charged reacting cases use zero. This does not certify physical convergence.
Source archives record backend, precision and gmpy2/MPFR/GMP versions. This is
an opt-in performance alternative, not a repair of physical convergence errors.

`interpolate` uses `precise_search.searchsorted`: monotone rounded values bound
the search, then original precise comparisons resolve every equal-float tie.
The returned index is exact; geometry, fractions and physical coefficients are
unchanged. There is no rounded-coordinate membership approximation.
