# Acoustic bound on extreme species mass coordinates

The full all-heavy interface trajectory with 5 ns source intervals completed
10 ns, then stopped after 73 accepted internal steps in interval three. Its
error was `2 eigenvectors failed to converge`, not a negative temperature or
physical-mass rejection. The preserved attempt is in
outputs/reactive_material_argon_all_density_jumps_5ns; it is not a 20 ns success.

Reconstructing the last source return from its saved moments/precise quadrature
and inherited knots, then replaying all 73 recorded steps, reproduces every
admissible bound and full energy exactly. The internal prefix lasts
6.92770080723507e-11 s. Ar2 alone repeats the failure at that state. The full
accepted state and its 53-by-53 mass/stiffness matrices are saved in
outputs/reaction_acoustic_eigenfailure.

The original bound uses the largest eigenvalue of

    K v = lambda M v,
    K = sum_i [gamma (p_i+Q_i) h_i/H_owner(i)^2] b_i b_i^T,
    dt_acoustic = CFL / max(sqrt(lambda_max), deformation_rate).

M is the actual consistent physical relative-coordinate mass, b_i the nodal
incidence vector, h_i the fine width and H the actual velocity-element width.
Pressure/numerical stress, every coefficient and all physical masses remain.
K is positive semidefinite for nonnegative pressure/stress; mass positivity is
not established by a successful eigenvalue routine alone.

The actual mass diagonal spans 7.37550e100 and its measured condition number is
1.03369e101. Define D=sqrt(diag M), y=Dv and solve

    (D^-1 K D^-1) y = lambda (D^-1 M D^-1) y.

This invertible congruence leaves eigenvalues and physical dynamics unchanged.
The balanced mass condition is 3.00528. Normalize balanced stiffness by its
largest absolute entry for the float solve, restoring that scalar in lambda.
There is no regularization, removed mass/node, eigenvalue floor or driver retry.
Exactly zero stiffness contributes zero acoustic frequency. Other invalid or
unrepresentable problems still reject; no acoustic step is enlarged by policy.

| Method | Largest lambda (s^-2) | Relative difference from reference |
| --- | ---: | ---: |
| Original generalized gvx | Fails to converge | — |
| Balanced gvx | 1.4421428904532682e17 | 1.04814e-16 |
| Balanced gvd | 1.4421428904532682e17 | 1.04814e-16 |
| Independent 100-digit solution | 1.442142890453268311157027...e17 | Reference |

The independent calculation interprets the saved float matrices exactly in
100-digit arithmetic, balances them, performs high-precision Cholesky, and
solves the resulting full symmetric spectrum with mpmath 1.3.0. A recovered
physical eigenvector gives componentwise original-equation backward error
3.19425e-70. Scripts, companions and spectrum.json preserve this evidence.
See [LAPACK DSYGVX](https://www.netlib.org/lapack/explore-html/db/de3/group__hegvx_ga73a91a9940b0a1572748844e48401bf8.html)
for its solver/error contract; failure wording alone is not a physics diagnosis.

The precise solver now uses `mass_scaling.largest_mass_eigenvalue` with the
explicit balanced gvx calculation. The previous banded acceleration/return
scaling is unchanged. A prior separate within-call width-reuse change preserved
all 202 arrays in a complete matched 20 ns case, including every actual dt and
budget; that equivalence does not verify the new acoustic arithmetic by itself.

Fresh validation uses the original full 20 ns all-heavy/5 ns-source trajectory
and all six independent neutral-decay/charged-affine cases at 16/32/64 cells.
They run from frozen /private/tmp/frc-balanced-acoustic-20261017; the coupled
source SHA is e8aeb92c292f20d6626ba14eecd785edabe382cd68dd2d4768bddfc16adfde81.
At this checkpoint the runs are live, and only decay_n16 is complete. Its full
6.083704 us errors include front 1.87926e-13 m, velocity 1.83809e-9 m/s,
energy 1.28419e-13 and population 9.024e-15 relative. Do not label the remaining
cases or coupled run successful until terminal evidence exists. Current,
source-time, thermal and construction-level physical validation remain open.
