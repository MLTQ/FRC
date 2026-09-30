# Physical regular-fluid and current coupling

This component connects the verified regular-coordinate fluid representation to
the compatible magnetic/current geometry. It is an instantaneous force/induction
component, not yet the complete moving two-species equations or a formation
integrator. Both previous solver closures remain unchanged.

## Common physical fields

Use the actual state and kinetic energy in REGULAR_COORDINATE_NUMERICS:
U=(n,q=n*u,l=n*Omega,mz,Ce,Ci), x=r²/2, with no independent axis q/l.
Let S sample the element degree-p scalar polynomial and the real periodic
trigonometric interpolant. Let V sample the transverse velocity polynomial,
degree p−1 in the first element and degree p elsewhere. At common physical
quadrature with weights W,

    nh=S n, vr=r V(q/n), vtheta=r V(l/n), vz=S(mz)/nh.

The even axial Nyquist coefficient represents a cosine. Radial field cells
subdivide fluid elements, so their quadrature does not cross discontinuous
element boundaries. Density is checked at points and common quadrature;
positivity between every possible point is not proved. No cell-average density
conversion, density floor, velocity filtering or alternate continuity is used.

## Force and torque adjoints

Let Vi be the fluid nodal volume. For supplied physical force density f,

    qdot_i = [V.T W(r fr)]_i / (M Vi ri²),
    ldot_i = [V.T W(r ftheta)]_i / (M Vi ri²),
    mzdot_i = ni [S.T W(fz/nh)]_i / (M Vi).

The first two formulas act only on active positive-radius nodes. Axis slots
remain zero by definition; there is no radius floor. With the actual grid
energy gradients M ri² q/n, M ri² l/n and M mz/n, this gives

    sum Vi H_U dot Udot_force = integral vh dot f dV.

Partition of unity gives exact integrated axial force and angular torque:
sum M Vi mzdot=integral fz and sum M Vi ri² ldot=integral r ftheta. These are
actual momentum statements, not work-only corrections. Pointwise force
accuracy still needs independent convergence evidence. Variable-density
continuous kinetic energy is not identical to grid energy; the preceding
fluid study measures that distinction rather than silently changing H.

## Current inertia on this density

With P the compatible current/electric sampling map, use the physical reduced
mass mu=mi*me/(mi+me), electric charge e, and

    I=P.T W [mu/(e² nh)] P,
    Kj=integral mu |Jh|²/(2e² nh),
    grad_n Kj=S.T W[-mu |Jh|²/(2e² nh²)].

The supplied n_t is the actual regular-fluid continuity derivative. Its
sample is S n_t, and its Ohm forcing is
F_n=−P.T W[mu Jh (S n_t)/(e² nh²)]. The same physical density map thus appears
in force, inertia, energy and changing-density work. Physical electron mass
is retained. The old original-equation current, Ampere and induction routines
are reused with this explicitly different density representation.

The complete affine field B=B0 ez+deltaB supplies Lorentz force. Ampere acts
on deltaB without loss of its small coefficients; local magnetic energy keeps
the bias cross term. The fixed-flux periodic geometry and natural Btheta wall
condition are unchanged. No vacuum domain or external coil boundary is added.

## Ideal/Hall pairing and the explicitly incomplete remainder

The implemented electromagnetic force is f=Jh×Bh. Electric forcing is
F_ideal=P.T W(−vh×Bh), F_Hall=P.T W[delta Jh×Bh/(e nh)], plus F_n, where
delta=(mi−me)/(mi+me). Hall work vanishes and ideal work exchanges exactly
with grid bulk kinetic work. The original finite-inertia solve still checks
Ohm, Ampere and Faraday residuals, including the existing explicit curl-kernel
conditioning diagnostics. No tolerance or physical mass is changed.

For this deliberately incomplete pair, the total energy rate is

    d(Hfluid + HB + Kj)/dt = integral mu |Jh|² (S n_t)/(2e² nh²).

That nonzero term is returned as `unresolved_density_exchange`. It is not
external heat, a numerical energy repair or evidence of total conservation.
The manufactured induction test has additional declared forcing; the ideal/
Hall check is a separate instantaneous evaluation without that forcing.

## Next physical transport contract

The continuum equations require the drift stress T=mu*n*w*w, w=J/(en), in
bulk momentum, and mu/e times
(v·grad)w+(w·grad)v−delta(w·grad)w in Ohm's law. Cylindrical vector derivatives
include the geometric terms. Using n_t=−div(nv) and div(nw)=0, their energy
identities are

    −mu*n*w·(v·grad)w + (mu/2)*w²*n_t = −div(Kj*v),
    v·(−div T) = −div(T·v) + T:grad v,
    −mu*n*w·(w·grad)v = −T:grad v,
    delta*mu*n*w·(w·grad)w = delta*div(Kj*w).

Thus transport, stress work and density exchange belong to one common discrete
derivation. A fitted force that merely cancels the measured global remainder
would not verify these equations. Next implement their actual shared physical
derivatives/adjoints and full species entropy/enthalpy/Poynting transport, then
independent original parent rates and unforced wave/contact evolutions.
Heated moving exact vacuum, physical closure/ionization, real boundaries and
converged FRC formation remain open. No numerical diffusion is added here.
