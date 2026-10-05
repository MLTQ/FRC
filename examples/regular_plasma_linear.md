# Separate stationary parent linearization

Uniform n, uniform individual entropies, zero flow/current and constant axial B
give independent linear equations: n_t=-div(n0 v), M n0 v_t=-grad(dp)+J cross B0,
C_s,t=s_s n_t and finite-mass ideal/Hall induction. Here
dp_s=(gamma_s-1)h_s[dC_s+(gamma_s-s_s)dn]. Current advection, drift stress,
changing-density inertia and current entropy are quadratic about this state.
Pressure enters E through sum a_s grad(dp_s)/(e n0), an exact compatible curl
kernel at linear order; it is retained in full E.

`StationaryLinearPlasma` directly implements centered regular-coordinate
radial/axial linear fluxes, dependent axis trace, shared pressure interfaces,
zero wall mass flux and original ideal/Hall plus physical current inertia.
It never differences the nonlinear RHS to construct its operator. Nonlinear
centered differences of mixed perturbations must independently converge to it.

This avoids feeding nonlinear O(epsilon^2) nearly curl-free pressure remainders
into an otherwise zero induction problem. It does not repair or change the
nonlinear solver's documented strict rejection of some isolated perturbations.
The class is a verification oracle at uniform rest, not a nonlinear time method.
