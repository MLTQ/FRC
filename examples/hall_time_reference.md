# Semidiscrete Hall time reference

Independent helicity decomposition and a dense matrix exponential solve the
spatially discretized linear system exactly to numerical exponential precision.
The two-by-two equations are assembled directly from the physical coefficients,
radial FE eigenvalues and axial wave numbers, without using the production
four-by-four operator or midpoint update. Recombining both helicities preserves
the actual interpolated initial data, including unintended secondary modes.

This is a temporal oracle for fixed spatial matrices, not a continuum reference.
It exposes slow-branch fast-mode contamination which an exact single-mode time
test does not contain. The Bessel reference remains the independent continuum
spatial and dispersion check.
