# Derivatives of conservative regular moment fields

Differentiates the same two mirrored degree-3/5 cell-containing reconstructions
as `compatible_velocity`. Radial fits use actual annular means of r^p P(r²/2)
for p=0,1,2. The target derivative includes both the r^p prefactor and the
polynomial derivative. Axial fits use physical periodic cell means.

These are interior physical derivatives of piecewise polynomials. They do not
include distributional jumps at cell boundaries and do not replace the shared
conservative continuity divergence. No derivative is inferred from fitted
current, analytic velocity or a density floor. Polynomial and independent
smooth-vector references must check axis and curl accuracy.
