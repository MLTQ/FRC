# Exponential local field/current evolution

Applies exp(dt A) to the same frozen-population skew Lorentz/Ampere system
documented in atomic_electromagnetic.md. Species momentum coordinates and electric
coordinates have the same quadratic-energy units. Diagonal species blocks map
h to h cross Omega_s; off-diagonal blocks are +omega_s I and -omega_s I.
Thus A.T=-A. The continuous finite-time map preserves the sum of their squared
norms. Its sparse representation contains the real physical frequencies, with
no artificial mass, phase fitting, damping or conservation projection.

The caller passes every occupied charged species, including arbitrarily small
nonzero populations. Exactly empty and uncharged populations are uncoupled and
remain outside this solve. This is not a threshold or omitted trace species.
Sparse matrix-vector products retain the small coupling factors in their own
rows, avoiding the contamination possible from a dense eigenvector basis for
extremely disparate populations.

Uses SciPy's
[expm_multiply](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.expm_multiply.html)
action and the analytic zero trace. The matrix exponential is accurate to the
library's floating-point algorithm, not symbolic/exact arithmetic. Large
frequency times dt can be expensive. The caller measures energy and phase
against independent analytic modes; neither is overwritten to pass.

This removes implicit-midpoint phase dispersion within a constant-population
source substep. Reaction/collision splitting, changing density and spatial
Maxwell coupling can still be underresolved. It does not establish that a large
global timestep resolves the complete plasma, and does not replace the moving
solver's quasineutral field constraint without a consistent new assembly.
