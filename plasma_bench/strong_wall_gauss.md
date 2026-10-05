# Gauss time evolution in the admissible wall space

Reuses the complete frozen ImplicitPlasmaGauss.step nonlinear algorithm,
strict stage/endpoint residuals, positivity checks, Newton fallback and atomic
failure. Only the constructor supplies StrongWallCoordinates and the restricted
fixed-background approximate inverse. Every complete RHS evaluation uses the
new explicitly declared wall Galerkin equations with unchanged physical masses.

The same fourth-order two-stage coefficients and tolerances apply, now in the
documented reduced physical basis. No time retry, artificial damping, density
floor or accepted-state projection is introduced. Boundary reactions and the
original unprojected curl defect require separate physical verification.
