# Physical quadratic perturbation-energy scale

For the physical conservative thermal variables (n,C=n*s), internal energy is
U=A*n^gamma*exp(C/n). At a positive background let a=delta n/n and
b=delta s=(delta C-s*delta n)/n. Half the Hessian contraction is

    U2 = U/2 * [(b+(gamma-1)*a)^2 + (gamma-1)*a^2].

The completed-square form stays positive for gamma>1 and avoids cancellation.
The independent high-precision test differentiates the original exponential
EOS along conservative n,C; it does not differentiate this expression.

`densities` adds the original stationary-background bulk kinetic, separate
thermal, physical current-inertia and magnetic perturbation energies. The full
harmonic bias has a linear cross term and no additional quadratic term. All
physical masses and fields are supplied explicitly. This is a normalization
scale, not a new evolution energy or a state/heat correction.

A growing mode on a nonuniform background need not conserve this positive
quadratic energy separately: nonlinear mean changes exchange background energy.
Only actual total-energy drift is divided by the initial quadratic scale.
No stability, confinement, free-energy invariant or formation claim is implied.

Nonfinite density, physical perturbations or masses reject; the metric never
turns an invalid physical input into a finite reported normalization.
