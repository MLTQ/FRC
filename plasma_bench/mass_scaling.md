# Physical mass-system scaling

`solve_mass_banded` solves the existing symmetric tridiagonal equation
`M v = b` after the invertible change of coordinates
`D = sqrt(diag(M))`, `A = D^-1 M D^-1`, `y = D v`.
The balanced equation is `A y = D^-1 b`; the returned velocity is `D^-1 y`.
This is the same equation, mass distribution and velocity space. It adds no
mass floor, diagonal regularization, population cutoff, viscosity or coarsening.

The actual rejected argon reaction projections had mass condition numbers of
about 10^29 and 10^46, while diagonal equilibration gave approximately 9.39
and 3.05. The unscaled band solve returned enormous velocities and negative
recovered thermal energy. Scaling and a separate 100-digit tridiagonal solve
both permit those same source snapshots to project without thermal repairs.
See `outputs/reactive_material_argon_subcycled_{10ns,5ns}` for the original
rejections and arithmetic diagnostics. Snapshot success is not a completed
charged evolution or a source-interval convergence result.

Inputs are finite SciPy tridiagonal band storage and a vector or matrix load.
Every diagonal entry must be strictly positive. Nonfinite/unrepresentable
results reject; the helper does not hide singular physical matrices. Sequential
division avoids unnecessarily underflowing products of diagonal scale factors.
Physical fields and coefficients remain float64.
