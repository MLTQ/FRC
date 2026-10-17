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

The precise reaction-return and material-acceleration paths use the band helper.
Their mass matrices and right-hand sides keep the original equations; source
projection and subsequent motion must both avoid the demonstrated unscaled
solve failure. Other historical solvers are unchanged.

`largest_mass_eigenvalue(K, M)` now applies the same congruence to the precise
acoustic generalized eigenproblem. With D=sqrt(diag M), solve
(D^-1 K D^-1)y=lambda(D^-1 M D^-1)y. Normalize the transformed stiffness by
its maximum absolute entry and restore that positive scalar in lambda. The
explicit `gvx` driver selects the largest eigenvalue. There is no driver retry,
matrix regularization, eigenvalue floor, removed species or changed time bound.
An exactly zero stiffness gives zero acoustic frequency; nonfinite transformed
matrices or a nonpositive/unrepresentable largest eigenvalue reject otherwise.
Positive definiteness of the full mass remains the eigensolver's requirement.
Inputs are the symmetric physical M and positive-semidefinite K assembled by
the compatible material equations; symmetry is a caller contract.

The actual full all-heavy 5 ns-source trajectory stopped after 73 accepted
steps in its third interval. Exact replay of all saved bounds/energies isolates
Ar2's 53-by-53 problem. Its mass condition is 1.03369e101, versus 3.00528 after
balancing. Original `gvx` repeats its convergence failure; balanced `gvx` and
`gvd` both give lambda=1.4421428904532682e17 s^-2. Independent 100-digit Cholesky
and symmetric eigensolution give 1.442142890453268311157...e17, a relative
difference of 1.05e-16 and physical componentwise backward error 3.19e-70.
See outputs/reaction_acoustic_eigenfailure. This arithmetic evidence does not
replace the required corrected full 20 ns trajectory or local convergence.
