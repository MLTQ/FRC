# Separate fluid refinement of the same extracted fields

Keeps the original field mesh and magnetic direction exactly fixed. Embeds the
original conservative fluid tangent by exact element subdivision to4/8/16
elements, preserving one-sided traces and checking physical tangent agreement
at independent quadrature16/20. Radially constant stationary base density is
required. This is a separate spatial diagnostic, not an evolved trajectory.

Evaluates both signs of the original full RHS at epsilon1e-4/1e-5 on each mesh.
Compares physical projected growth and nonparallel rate with the original
direction and the original saved RHS. Saves all coefficients, solver residuals,
physical integration controls, hashes and source; fresh replay is required.
Refining one representation cannot establish full spatial convergence, and
unchanged jumps in the original direction remain part of the diagnostic.
