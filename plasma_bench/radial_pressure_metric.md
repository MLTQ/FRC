# Diagnostic physical radial pressure metric

Given bulk thermal plus kinetic face-work coefficients w_r, pull them back through the exact
normal particle-momentum map B_r. Solve

    [R_r.T W (M/n) R_r] m_r,t = −B_r.T w_r.

The matrix is the true physical radial kinetic Hessian at fixed actual n.
Its contraction with m_r therefore cancels the supplied radial face power.
Preconditioning uses the separable constant-density physical Gram matrix;
normalization preserves the actual total species mass and variable density.
Failed relative residual or nonfinite solve rejects without modifying state.

This is a diagnostic radial force, not a conservative Maxwell/drift stress or
a full energy theorem. It does not alter axial/angular momentum, and cannot
fix the separate uniform-translation obstruction of the axial thermal update.
Its independent radial force/equilibrium accuracy must be verified before
adoption. It is not called by the frozen baseline LogFluid study.
