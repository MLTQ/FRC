# Velocity accuracy with conservative energy constraints

`constrained_velocity` is an explicit research alternative for reaction return.
It minimizes the occupied-volume quadrature norm of v_h-u_source rather than
the density-weighted norm. This gives dilute material equal importance per
occupied volume, preventing a dense region from determining a huge extrapolated
velocity in a trace tail. It retains the same P1 knots and P0 amounts, without
mass floors, deleted species, coarsening or modified physical source moments.

The optimization has the exact physical constraints

    sum_i a_i v_i = P_source,  a_i = integral rho_target phi_i dx,
    K_c(v) = (1/2) integral_c rho_target |v_h|² dx <= E_source,c.

Thus momentum is conserved per connected species block, and U_c=E_c-K_c
preserves each fine cell's total energy. The constraints are convex quadratics;
the volume objective is quadratic. Heat is recovered once from the solved
velocity; there is no clipping, assigned floor or energy reset. Local momentum
redistribution and thermal changes are numerical projection effects.

First solve the equality-constrained volume problem analytically with two
banded solves. If its cell energies are admissible, no nonlinear optimization
is performed. Otherwise the current implementation requires a strictly feasible
constant-mean velocity start. If unavailable, it rejects rather than inventing
a feasible state. Exactly cold, zero-energy donors are not generally supported
by this interior-start strategy; the default solver remains unchanged.

The logarithmic barrier is -mu/ncells sum log((E-K)/E), with normalized velocity
coordinates and mu=1e-4,1e-6,1e-8,1e-10,1e-12. These are optimization tolerances,
not physical temperatures. Each Newton direction preserves the momentum equality;
line search remains strictly feasible and uses Armijo decrease. Each barrier
stage permits 80 iterations; a failure is terminal, with no timestep retry.
The maximum normalized momentum residual must be <=1e-12. The remaining barrier
parameter and iteration counts are reported. This procedure approximates the
constrained optimum from its interior; barrier sensitivity remains to measure.

The objective uses original occupied source quadrature points; no velocity is
invented at empty source points. A deficient quadrature matrix rejects. A
representable admissible velocity is unchanged by the equality solve up to
roundoff in the constant-density, zero-source limit. Nonlinear reacting states,
current/front accuracy, full time/space refinement and physical validation remain
required. This is not a general invariant-domain transport theorem.

Contract: receives precise model geometry, flattened source points/weights,
relative source velocities/occupancy, and normalized physical block momentum
and cell energies. Returns the same nodal velocity shape and optimization
diagnostics. Caller must recover and validate actual thermal energy using its
normal physical-energy expression, not the optimizer's normalized slack.

With an explicitly extended material basis, source quadrature includes the entire physical support and clamps element ownership at both edges. It integrates affine extrapolation shapes there, preserving the same momentum and fine-energy constraints.
