# Finite-pressure cylindrical Hall reference

This reference supplies an independent oracle for the next nonlinear solver.
It is a fully filled, singly charged, cold-electron (pe=0), massless-electron
Hall fluid with adiabatic ion pressure, constant background rho/p/Bz, zero
resistivity, periodic z and a conducting slip radial wall. It is not a
plasma/vacuum interface or a validated argon formation closure. Its waves are
linear; their finite-amplitude insertion has quadratic nonlinear residuals.

The species starting point is ion momentum, electron force balance, Ampere and
quasineutrality, with vi=v and ve=v-J/(ne). The massless limit of the species
equations removes electron drift energy and the extra drift momentum stress;
these cannot be omitted in a finite-mass extension. See the parent equations
in [Andres et al., sections II–III](https://arxiv.org/pdf/1410.6518).
The compressible adiabatic pressure law is an additional declared closure here.

## Continuum reduction

Let qR be a positive J1 zero, k=2 pi ell/L, K²=q²+k², and use exp(ikz).
The complex amplitudes are

    Br=-ik a J1(qr), Bz=q a J0(qr), Btheta=t J1(qr)
    vr=u J1(qr), vtheta=w J1(qr), vz=v J0(qr), p1=p J0(qr).

The equations derived directly from momentum, pressure and Faraday are

    rho u_dot = q p + B0 K² a/mu0
    rho w_dot = ik B0 t/mu0
    rho v_dot = -ik p
    p_dot = -rho cs² (q u+ik v)
    a_dot = -B0 u - ik DH t
    t_dot = ik B0 w - ik DH K² a,

with cs²=gamma p0/rho and DH=B0/(mu0 ne e). Here DH has units m²/s and
keeps the bias sign. rho1=p1/cs²; the independent stationary entropy mode is
excluded from these six propagating amplitudes.

In coordinates y=(u,w,v,p/(rho cs),K a/sqrt(mu0 rho),t/sqrt(mu0 rho)), the
generator A is anti-Hermitian. Thus i A is Hermitian and its eigenvalues are
real signed frequencies omega for exp(ikz-i omega t). Branch indices order
the three positive frequencies; they do not label a fixed curl helicity.
Both axial signs, both frequency signs and signed bias are retained.

Independent elimination gives

    (omega²-k² vA²)
      [omega⁴-K²(cs²+vA²)omega²+k² K² vA² cs²]
    = omega² DH² k² K² (omega²-K² cs²).

At DH=0 this gives Alfven and fast/slow magnetosonic waves. At q=0 it gives a
Cartesian sound wave and the two parallel Hall branches. q=0 is not a
nontrivial cylindrical J1 mode. In the cs→infinity limit the acoustic pair
separates and the remaining branches recover the earlier incompressible
quadratic. Finite cs does not justify using that older quadratic unchanged.

## Independent tests and energy

At a J1 root, integrals of r J0² and r J1² agree: I=R² J0(qR)²/2. The real
periodic wave has quadratic perturbation energy

    W = (pi L I rho / 2) |y|².

This includes p1²/(2 rho cs²), kinetic and magnetic perturbation energies.
It is distinct from the background total energy and is checked by physical
volume quadrature at different phases. All three branch fields satisfy
independently differenced Ampere/Faraday, electron force balance, ion momentum
and pressure equations. The dispersion polynomial and limiting formulas are
separate checks, not a copy of a numerical Hall evolution matrix.

Br, Btheta, vr and vtheta vanish at the wall. Ideal Hall tangential E vanishes
there, while the axis is regular. Jz generally does not vanish at the wall:
eta Jz would violate a resistive conducting boundary. These modes therefore
do not verify resistive damping. Tests use a separately stated pe=0 reduction;
warm barotropic electrons require their energy/enthalpy treatment even when
their pressure gradient makes no induction curl.

`examples/compressible_hall_reference.py` contains this oracle.
`hall_work_snapshots.wave_hall_rate` separately differentiates only the Hall
part of induction, not the full wave time derivative. The 42-case spatial
study uses it for current/Faraday consistency. Full nonlinear wave evolution
and time-step convergence have not yet been implemented for the new operator.
