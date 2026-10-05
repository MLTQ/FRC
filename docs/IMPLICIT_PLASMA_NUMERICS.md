# Full-equation implicit plasma time integration

This is a bounded verification extension of RefinedRegularPlasma. It advances
the original finite physical-mass n,q,l,mz,Ce,Ci and magnetic equations together
on a positively filled periodic cylinder. The prescribed harmonic bias is fixed.
It adds neither exact vacuum nor a laboratory formation model. Original spatial
fluxes, species/current work, physical constants and strict induction checks
are unchanged. Older source closures remain frozen.

## Time equations and acceptance

Implicit midpoint solves delta=dt F(U0+delta/2) and returns U0+delta. Two-stage
Gauss collocation solves Yi=U0+dt sum_j Aij F(Yj), then returns
U0+dt(F(Y1)+F(Y2))/2, with A=[[1/4,1/4-sqrt(3)/6],
[1/4+sqrt(3)/6,1/4]]. F is the complete actual nonlinear plasma RHS at each
stage, including density-dependent current inertia and both species laws.

Both solve in fixed physical active coordinates. The initial volume means
define n0, separate pressures and c²=sum gamma_s p_s/(M n0), M=mi+me.
Scales are n0 for n/Cs, n0 c/R for q/l, n0 c for mz,
sqrt(mu0 M n0)c R² for psi and sqrt(mu0 M n0)c/R for g. Only the already
constrained axis q/l and axis/wall magnetic flux slots are removed. Scales
remain fixed across a trajectory; a solver initialized at a later checkpoint
does not necessarily have the original coordinate/preconditioner contract.

The original residual must satisfy the active-coordinate infinity norm
atol+rtol max(|increment|inf,|dt times stage RHS|inf). Default atol=1e-12,
rtol=1e-10; wave accuracy controls explicitly tighten these. Gauss checks actual
rounded stage arrays and endpoint rounding. Midpoint reconstructs and rechecks
the actual returned rounded increment. Every RHS also retains every unchanged
original induction residual bound, independently of this time-equation norm.
Both point and common-quadrature states require positive finite physical EOS.
Failure leaves the caller's accepted state and magnetic arrays unchanged.

## Iteration and preconditioning

An explicitly approximate uniform physical background provides a sparse
regular-fluid matrix, a compatible physical-inertia/Hall magnetic block and
reciprocal ideal-induction/Lorentz work blocks. Uniform-flow fluid linearization
is separately derived and checked by nonlinear differences. At rest the full
matrix is checked against the existing independent parent linearization. The
uniform background never replaces a stage state or evolution equation.

A real or complex shifted block Schur inverse gives a chord correction.
For Gauss, a conjugate pair reduces the real two-stage inverse to one complex
shift h=dt(1/4+i sqrt(3)/12). The actual plasma state stays real. A line search
requires decreasing original residual or satisfaction of its unchanged bound.
Invalid physical-state trials are damped, up to 13 halvings; strict inner
induction failures propagate. No time retry is hidden in the stepper.

When that chord direction cannot decrease the residual, centered products of
the same full residual give a matrix-free Newton/GMRES correction. The physical
coordinate perturbation is cbrt(machine epsilon)(1+|increment|inf). GMRES uses
the same approximate inverse, relative linear tolerance 1e-4 and at most five
cycles of 30 directions. The proposal must still pass the original line search,
positive-state checks and final nonlinear bound. Linear iterations, fallbacks,
backtracks and all original induction maxima are reported. Changing the method
of solving the equation is distinct from changing an acceptance tolerance.

Dense fluid/field and magnetic blocks limit this implementation to bounded
verification meshes. A scalable formation solver will need a structured or
matrix-free block preconditioner. Performance on this reference is not a claim
of usable large-device throughput.

## Accuracy, diffusion and limits

Midpoint is second order, R(z)=(1+z/2)/(1-z/2). Gauss is fourth order,
R(z)=(1+z/2+z²/12)/(1-z/2+z²/12). Both have unit modulus on the imaginary
axis and are A stable, not L stable. They supply no intentional damping to
linear oscillatory modes; a large step still gives wrong unresolved phase.
These properties do not prove stability for time-dependent nonlinear states.
No viscosity, filter, artificial mass, numerical atmosphere, temperature floor
or total-energy repair is introduced. Nonlinear physical energy is not exactly
quadratic and is not exactly conserved by either time method.

Complete contact/wave trajectories require independent time/mesh refinement,
full E, both temperatures, B/J, point/sampled positivity, axis/wall/max errors
and grid versus reconstructed continuous energy. Strong growth in a known
contact must remain visible even if all nonlinear equation residuals pass.
The captured chord-only failures and their accepted checkpoints are retained;
a corrected solver iteration alone does not validate long-time physics.

## Radau control and spatial diagnosis

A separate two-stage Radau IIA implementation uses A=[[5/12,-1/12],
[3/4,1/4]], b=[3/4,1/4]. The second physical stage is the returned state,
and its original weighted endpoint equation is independently checked. The
complex shift is dt(1/3+i sqrt(2)/6); all strict nonlinear/induction checks and
Newton fallback are shared in meaning with the preceding methods.
R(z)=(1+z/3)/(1-2z/3+z²/6) is third order and L stable. It deliberately
damps unresolved oscillations, including physical fast waves; that temporal
dissipation is measured and must not be called physical resistivity.

Radau also develops long-contact growth. A stationary-contact spatial Jacobian
now exposes rapid growing modes independently of time integration. A separate
original-equation Fourier/Bessel linearization has slower growing modes in its
checked subspaces. Those radial spaces are not yet matched. Read
IMPLICIT_PLASMA_VERIFICATION.md and stationary_contact_continuum.md before
interpreting either result; time damping has not solved this coupling gate.
Ranocha and Gassner's Euler density-wave study motivates checking local linear
stability independently of conservation, but does not diagnose this plasma:
https://arxiv.org/abs/2009.13139 .
