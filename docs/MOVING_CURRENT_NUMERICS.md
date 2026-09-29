# Moving finite-mass induction and conservative continuum contract

`MovingCurrent` evaluates an instantaneous induction rate for prescribed positive
density, bulk velocity, density derivative and species thermal state. It retains
physical electron/ion mass. It does **not** yet evolve those quantities together
or replace the existing guarded formation path. The verification domain is a
filled periodic cylinder; it has no physical plasma–vacuum edge.

## Parent species and exact moments

Use the SI conventions of CURRENT_THERMAL_NUMERICS:

    M=mi+me, rho=M n, alpha=me/M, beta=mi/M,
    delta=beta-alpha, mu=mi me/M, w=J/(ne),
    vi=v+alpha w, ve=v-beta w, J=curl B/mu0.

Both species obey continuity, isotropic pressure/internal energy and momentum:

    n_t+div(n vs)=0,
    ms n (partial_t+vs dot grad)vs = qs n(E+vs cross B)-grad ps+Fs,
    (ps_t+vs dot grad ps+gamma_s ps div vs)/(gamma_s-1)=Qs.

Production forcing is absent unless explicitly supplied. The independent
manufactured reference computes Fs and Qs from these original equations, with
qs=+e/-e. Its nonzero sources are part of the declared reference, not a fitted
discrete residual or an unforced equilibrium.

Quasineutrality, no displacement current and div J=0 give a common continuity
equation. Summing and differencing species momentum gives

    rho (partial_t+v dot grad)v
      = J cross B-grad(pi+pe)-div(mu n w w)+Fi+Fe,

    E = -v cross B+delta w cross B+(alpha grad pi-beta grad pe)/(ne)
      +mu/e [w_t+(v dot grad)w+(w dot grad)v-delta(w dot grad)w]
      +(beta Fe-alpha Fi)/(ne).

For zero external force the conservative momentum flux tensor is

    rho v v + mu n w w + (pi+pe+B²/(2mu0)) Identity - B B/mu0.

In particular, the bulk inertia is **(mi+me)n**, not mi n. Drift stress belongs
in this flux and its cylindrical geometric terms. It must not be appended as
an unrelated heating correction. Pressure anisotropy, kinetic stress, charge
separation, collisions, ionization and conduction are not included here.

The need for current kinetic energy and the corresponding momentum term is
discussed by [Kimura and Morrison (2014)](https://arxiv.org/abs/1406.2745).
Our exact mass fractions, moving-density transformation and numerical method
are derived here and checked against both parent species; the paper's reduced
order approximations are not silently substituted for those equations.

## Changing density and velocity gradients in induction

Let B*=B+(mu/e)curl w. The identities

    (v dot grad)w+(w dot grad)v
       = grad(v dot w)-v cross curl w-w cross curl v,
    (w dot grad)w = grad(|w|²/2)-w cross curl w,
    w_t=J_t/(ne)-w n_t/n

produce the equivalent equation (with the force source added when present):

    E = (delta w-v) cross B* - (mu/e) w cross curl v
        - (mu/e) w n_t/n + mu/(ne²) J_t
        +(alpha grad pi-beta grad pe)/(ne)
        +grad[(mu/e)(v dot w-delta |w|²/2)].

The supplied n_t must ultimately come from the same conservative particle flux
as the bulk continuity update. Merely rebuilding an inertia matrix at a changed
density omits this term and the velocity-gradient term. The present manufactured
study supplies the analytic n_t, so it does not discharge that integration gate.

`dual_kinematics` interpolates the explicitly regular nodal components
(ur=vr/r, Omega=vtheta/r, vz) as Q1 functions in x=r²/2,z. At quadrature:

    curl v = (-r Omega_z, r(ur_z-vz_x), 2 Omega+r² Omega_x).

Axis limits follow from these regular variables without a division floor. A
conservative primal-momentum to this dual representation map remains to derive.
Current traces and the kinetic/mixed scalar potentials use the same stated
interpolation as the fixed-density component; measured interpolation error
remains and is not removed by the exact commuting gradient.

## Instantaneous linear solve and physical work

Write C=-curl, H for the compatible magnetic Hessian, Me for the electric mass
and I for the physical current-inertia matrix. With f the complete reduced weak
electric forcing, the mixed equations are

    Me E = I J_t+f,    Me J_t = -C.T H B_t,    B_t=C E.

Fixed axis/wall flux rates are eliminated. On the remaining magnetic space,

    [H+H C Me^-1 I Me^-1 C.T H] B_t = H C Me^-1 f

is positive definite. Diagonal scaling balances flux and toroidal units; an
exact factorization of the scaled H is only a preconditioner. Conjugate
gradients solves the full operator. Original Ohm/Ampere/Faraday and operator
residuals are checked, including Faraday in the magnetic energy norm. Up to
three linear iterative refinements correct the original reconstructed residual
without changing the tolerance or equations. Failure leaves inputs unchanged.
These linear corrections are not hidden physical time-step retries.

The returned full E restores both pressure and kinetic/mixed scalar gradients.
The tangent work identity includes the actual density derivative:

    d(WB+KJ)/dt = -j dot f - (1/2) integral mu |J|²/(ne²) n_t/n dV.

An independent physical r,z quadrature differentiates magnetic/current energy
along B+h B_t, J+h J_t, n+h n_t. Symmetric energy differences are evaluated
before quadrature using polarization/common denominators, then extrapolated to
fourth order at two intervals. This avoids subtracting nearly equal global
energy totals; no state, pressure or energy is repaired.

## Full local energy and remaining discrete coupling

Define Ks=ms n |vs|²/2, Kbulk=rho |v|²/2, KJ=mu n |w|²/2,
epsilon_s=ps/(gamma_s-1), Hs=epsilon_s+ps. The exact species sum is

    W=B²/(2mu0)+Kbulk+KJ+epsilon_i+epsilon_e,
    F=(Kbulk+KJ+Hi+He)v + E cross B/mu0
        +w[mu n v dot w-delta KJ+alpha Hi-beta He],
    W_t+div F = vi dot Fi+ve dot Fe+Qi+Qe.

Independent direct species kinetic/enthalpy fluxes test this reduction. The
mixed `mu n (v dot w) w` transport term is essential; an omitted-term control
retains a finite local-energy defect under refinement.

The next implementation must advance particle number, conservative bulk
momentum/angular momentum, both species thermal quantities and full total
energy with shared fluxes and consistent temporal means. It must retain both
component-energy and conservative-cell-energy accounts and measure their local
discrepancy; resetting one to match the other is not a conservative solution.
Verify finite-mass compressible waves and this moving manufactured trajectory,
spatial/time convergence, pressure positivity, local work and entropy behavior.
An instantaneous work identity does not prove a conservative finite time step.

## Boundaries, numerical error and scope

The cylinder has periodic z, fixed psi on the axis/wall, Btheta(R)=0 and free
Ez(R). Bulk radial flow is zero at the wall. These are bounded reference
conditions, not a conducting Ez=0 wall, real gas ends or vacuum matching.
No numerical atmosphere is introduced. The physical exterior is not modeled
by this component; eventual formation still requires exact vacuum outside.

No viscosity, resistivity, filter, upwind diffusion, density floor or artificial
mass is added to this rate. Finite element representation, projection, nodal
potential interpolation and finite solve tolerance produce the measured
spatial errors. There is no time integrator here, hence no temporal dissipation
or stability claim. Actual fluid integration, shocks/entropy, physical transport,
the heated free surface, exact vacuum, real ends/coils and kinetic validity
remain required before formation experiments can be trusted.
