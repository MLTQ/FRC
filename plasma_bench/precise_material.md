# Compatible material dynamics with explicit precise geometry

`PreciseMaterialSubcells` retains the existing fine P0 density/thermal and nested
P1 velocity equations, consistent mass, physical species scaling, minmod Q,
pressure/work pairing and electric quadrature. It inherits the existing fine
velocity, momentum-load restriction and factored number/energy conversions.
Only geometric arithmetic changes: supports, reference fractions and fine-face
velocities carry50-digit decimal coordinates; width/shape ratios are formed
before conversion to float physical coefficients. No density cut is removed.

`PreciseMaterialState.added` advances geometry within that context. The full RK4
composition in `PreciseMaterialSet.step` uses the same precision for stage sums.
Retaining fine-face velocity differences matters as well as accurate position
addition: rounding those velocities first can erase the represented gradient.
Velocity, thermal state, mass solves, field coefficients and budgets remain
float64. This is not an arbitrary-precision physical-state solution.

`PreciseMaterialSet` keeps the original distinct-support validation and explicit
acoustic/plasma restriction. Union sorting, interval membership and quadrature
use precise positions. Gauss E and its energy use float physical lengths after
exact-coordinate subtraction, including the empty exterior. Sampling and totals
retain the laboratory frame, actual currents and all scaled species inventories.
Negative heat, nonpositive/unrepresentable widths, overlapping same-species
blocks and contact with remote ends still reject. No repaired widths or energy.

This explicit alternative permits continuing precision work while the earlier
scalar-geometry runs keep their frozen source. It is not formation adoption.
The original equations/assumptions are in material_subcells.md and
material_species_set.md. Source projection and saved geometry must use the
corresponding precise path; casting these faces back to float loses the fix.
