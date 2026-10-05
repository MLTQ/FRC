# Faster unchanged bounded induction equations

UniformWallInduction changes only the inner iterative approximate inverse from
H to the fixed uniform-density physical current-inertia metric. It inherits
every original wall equation, analytic kernel separation, actual changing
density, residual, refinement limit and atomic rejection. No solution/rate is
fitted or corrected. All other plasma terms and normal-wall reaction are inherited.

UniformWallPlasma caches the density-independent inverse and operators, rebinding
their *actual* model for every RHS evaluation. `preconditioned(system)` makes a
new assembly on the same mapping, geometry, fluid, masses and tolerances; it does
not mutate the original system or any state. Cache-disabled and variable-density
controls must remain physically equivalent within the unchanged strict bound.
This remains a filled-periodic verification model, not exact vacuum.
