# Finite-cylinder vacuum response to an independent dipole

Solves the actual new magnetic-only exterior with the boundary trace of an
axially displaced point dipole located inside the excluded plasma cylinder.
The dipole is a mathematical field fixture, not a point current placed into a
physical plasma simulation. Exact flux and field are analytic Cartesian dipole
formulas expressed in cylindrical coordinates, independent of the production
multipole/electromagnetic sampling functions.

Compares physical field around the side, at both ends, on the external axis and
beyond the outer matching rectangle. Reference exterior energy is independently
integrated on the inner cylindrical boundary using the exact dipole field and
Green's identity. The solved exterior includes its all-space magnetic energy.
Patch and multipole errors are measured, not inferred from matrix symmetry.
The test is a static physical solve supporting finite-end integration; it does
not replace the required full finite-coil/plasma evolution. No unit tests or
new reporting/replay framework. Saved coefficients permit direct field reuse.
