# Full Gauss solve with modal preconditioning

Reuses the frozen complete ImplicitPlasmaGauss.step algorithm, physical changing
density, separate species and full StrongWallPlasma equations. Only admissible
coordinate layout and the fixed-background proposal inverse differ from the
earlier dense StrongWallGauss constructor. Stage/endpoint equations, positivity,
strict induction, Newton fallback, backtracking and atomic failure are unchanged.

The residual uses the documented separable physical coordinate infinity norm;
compare actual reconstructed fields and time equations when assessing the new
layout. FFT operates only on preconditioner corrections, never as a filter on
the evolved state. No artificial damping, altered mass or energy repair.
