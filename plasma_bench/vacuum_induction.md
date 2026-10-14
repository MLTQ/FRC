# Inertial induction with free boundary flux

Uses the existing strict SplitInverseWallCore equations and residual bounds,
with a new admissible lift including every non-axis psi node, including R.
Toroidal field remains in the existing exact wall-trace null space. The positive
magnetic metric includes the actual eliminated vacuum energy. Ampere, inertia
and Faraday therefore share the same current and boundary degrees of freedom.

The constant-density current-inertia matrix is only a preconditioner, factored
once. Every nonlinear RHS uses actual density. Kernel separation, inner physical
and coefficient residuals, full Ohm restoration and difference-first full-trace
Faraday checks retain existing tolerances and iteration limits. No state repair,
physical mass change, extra diffusion or acceptance relaxation is introduced.
The full radial boundary current is determined by weak Ampere, not a separately
evolved auxiliary mode or prescribed E/eta replacement.

Five complete coupled runs and measured spatial/time/model comparisons are in
docs/VACUUM_CYLINDER_RESULTS.md. This bounded periodic coupling is complete;
finite coils/ends, moving free surfaces and formation remain open.

FiniteInduction supplies finite-axis lift/kernel maps through separate_kernel;
the strict solve and all residual bounds are shared unchanged.
