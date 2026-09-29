# Conservative electron pressure and entropy component

This component solves a missing part of the warm plasma equations: adiabatic
electron transport by the current, coupled to the electron-pressure term in
Ohm's law. It is an isolated operator with fixed positive density. It does not
yet advance ions, bulk momentum, Hall J cross B, resistivity, collisions,
conduction or vacuum. It must not be presented as the completed two-temperature
Hall solver. The prior cold-electron Hall model is unchanged.

## Parent split equations and energy

For singly charged massless electrons, quasineutrality, isotropic pe and an
adiabatic ideal electron gas, ve=v-J/(n e). During this split v=0 and n is fixed;
div J=0 preserves electron particle continuity. Retain

    E_p = -grad pe/(n e),  Bdot = -curl E_p
    epsilon_e_dot + div(epsilon_e ve) + pe div ve = 0,
    epsilon_e = pe/(gamma_e-1),   J=curl B/mu0.

This is a term decomposition of the parent equations, not a claim that the
other terms are physically absent over a real formation trajectory. In an
axisymmetric scalar thermal state, Jtheta does not advect entropy. E_p has no
theta component and changes only Btheta; poloidal energy stays constant.

Let a=(n/nref)^(gamma_e-1), w=pe/(n a), U=n w, A=a/(gamma_e-1). Then

    Udot = div(w J/e),  epsilon_e=A U
    E_p = -grad(pe/n)/e - w grad(A)/e.

The continuum electron energy flux is (epsilon_e+pe)ve. Combined with pressure
Poynting transport, total energy obeys

    (epsilon_e+B²/(2mu0))dot
      + div[(epsilon_e+pe)ve + E_p cross B/mu0] = 0.

Thus the entropy formulation includes the required electron enthalpy drift and
pressure work. Adding only pressure Ohm to magnetic evolution would omit this
energy transfer. nref is a positive proxy scaling, not a changed physical density.
The code verifies its invariance explicitly.

The parent species/thermal framework is described in the
[NRL formulary](https://www.govinfo.gov/content/pkg/GOVPUB-D210-PURL-LPS118567/pdf/GOVPUB-D210-PURL-LPS118567.pdf),
printed pp.36–39. The separate-temperature, finite-mass and edge reductions
still need the checks in TWO_TEMPERATURE_PLAN.

## Actual conservative charge flux

Set x=r²/2. Primal cell g=Btheta/r has magnetic energy

    H_B = sum_cell V <r²> g²/(2mu0),  <r²>=x_left+x_right.

Define a current stream function at primal cell centers:

    K_i,j = 2pi (x_i+x_i+1) g_i,j / mu0.

Dual scalar nodes lie at (x_i,z_j); their control faces pass through adjacent
primal centers. Interior radial current is K_i,j-1-K_i,j; axial current at
node i is K_i,j-K_i-1,j. These are **integrated amperes**, not sampled densities.
Each dual control volume has exactly cancelling outward currents, including
the regular axis and periodic seam. Its volume is 2pi delta(x_dual) delta(z),
with actual half-support volumes at the axis and wall.

At the wall use K_wall=2pi R Btheta_wall/mu0. Its radial difference supplies
the actual normal wall current; its axial difference from the last interior K
supplies the boundary dual current. Omitting the physical trace invents a
return-current sheet. For Btheta=C r the correct current is Jz=2C/mu0 at every
dual node, including the wall. The omitted-trace control gives the wrong sign
there and the independently known sheet offset.

The scalar gradient has Er=r delta(phi)/delta(x), Ez=delta(phi)/delta(z), and
Faraday gives gdot=delta(Ez)/delta(x)-delta(Er/r)/delta(z). Its curl annihilates
the scalar gradient. Against any electric field the exact work identity is

    H_B_dot = -sum_edges I delta(length)*E + sum_wall K_wall delta(z) Ez_wall.

The final term is inward magnetic Poynting power. Restricting Ez_wall=0 and
dropping its unused degree of freedom recovers the previous Hall adjoint.
No mass projection is used to manufacture divergence-free charge transport.

## Thermodynamic spatial and temporal means

Use the same edge w in U transport and pressure electric field. For endpoints
wL,wR, choose the reciprocal logarithmic mean of reciprocals:

    L(a,b)=(b-a)/(log b-log a),  L(a,a)=a
    w_edge = wL wR/L(wL,wR).

Then w_edge(1/wL-1/wR)=log wR-log wL. Local charge conservation makes both
sum(V U) and the discrete physical entropy S*=sum(n V log(w/e)) constant.
S* equals (gamma_e-1)S_e/kB up to an irrelevant constant entropy gauge. This
is a lumped nodal discretization of physical entropy, not a new equation of
state or proof that every reconstructed continuum entropy integral is exact.

The same w_edge makes global electron energy exchange cancel magnetic work.
Constant w (adiabatic barotrope) produces zero induction. Constant temperature
gives w=T/a, so w_edge delta(A)=T delta(log n); its pressure field is also
curl-free. The two distinct zero-battery limits are explicitly verified.
Stable near-equality series, relative atanh and log1p evaluations avoid losing
precision by subtracting logarithms of tiny SI temperature values. Independent
high-precision Decimal values test these means across extreme ratios.

For time evolution use g_mid=(g_old+g_new)/2 and U_stage=L(U_old,U_new).
Solve the coupled implicit equations with both spatial updates evaluated there.
The exact finite difference of log U is delta(U)/L(U_old,U_new), so the same
current-flux cancellation preserves finite-step S*. Magnetic midpoint preserves
quadratic H_B while electron energy is linear in U. The method is symmetric
and second-order consistent. This is a thermodynamic discrete-gradient update,
not arithmetic midpoint for U and not a post-step energy/entropy repair.

The general use of thermodynamic means in conservative flux design is discussed
by [Chandrashekar](https://arxiv.org/abs/1209.4994). That paper's compressible
Euler flux is not this cylindrical electron discretization; the derivation and
independent verification above establish this component's identities.

## Conditioning, boundaries and limitations

For pressure induction the exact curl-free grad(pe/n) part is omitted. Each
remaining closed A circulation also subtracts a local constant w, another
exact gradient null identity. This reduces cancellation near the axis without
altering the physical curl or fitting an axis multiplier. Full E_p is returned
for physical work and independent curl comparisons.

Picard iteration checks the actual final equations. State scales include initial
values and predicted changes, so a zero initial magnetic field is not hidden by
thermal background. Failed trials preserve inputs. There is no monotonicity or
unconditional-positivity theorem; nonpositive or nonfinite states fail explicitly.
No artificial heat, entropy source, temperature floor or diffusivity is added.
Sharp gradients/shocks will require an entropy-stable dissipative extension.

The current closed thermal reference requires constant wall Btheta and wall
density. Normal particle current vanishes, and the generally nonzero pressure
Ez at the wall has zero **integrated** magnetic input. Local boundary Poynting
flux need not vanish. This fixed magnetic-trace boundary is **not** HallMHD's
conducting Ez=0 wall. Deleting the gradient part for induction does not make
the physical electric field a conducting-wall field. General wall particle
inflow, surface charge, sheath and reservoir thermodynamics remain unresolved.

Next couple the dual thermal state to conservative moving ion mass flux and
the primal total-energy state, with explicit boundary compatibility and positive
ion/electron partition. See TWO_TEMPERATURE_PLAN; the isolated component must
not be installed by resetting primal energy or interpolating density arbitrarily.
