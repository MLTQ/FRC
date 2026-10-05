# Actual signed nonlinear continuum-mode trajectory

Evolves actual conservative fluid and admissible magnetic coefficients from a
signed independent Fourier/Bessel perturbation of the stationary contact.
UniformWallPlasma changes only an inner inverse; RecycledModalGauss advances
the full original bounded equations and physical masses. Every nonlinear stage
and endpoint must meet the explicit unchanged tolerance and positivity rules.
The independent mode never supplies a later state, forcing or fitted growth.

The caller declares separate meshes, signed amplitude, duration, step count
and absolute tolerance. Every accepted endpoint and actual stage change is
saved with a hash, full iteration history and induction diagnostics. Summary
progress is updated after each saved step. Exceptions preserve the last accepted
state and write a terminal failure; success requires every requested step.
No automatic time retry, residual relaxation or field/energy repair is used.

The source/runtime/companions freeze before evolution and must remain unchanged.
`states.npz` has exact initial/final states for independent quadrature audits.
Complete replay must also compare all per-step files and non-runtime metrics.
This periodic filled-cylinder verification is not FRC formation or a lab recipe.
