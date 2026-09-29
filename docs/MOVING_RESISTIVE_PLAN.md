# Next component: resistive material induction and heat

This is a derived implementation/verification plan, **not an implemented model
or completed verification result**. The current moving reference is ideal.
First test induction on frozen geometry; only then combine with material motion
and recheck source work, splitting error, pressure and interface accuracy.

## Fixed-geometry compatible form

Let K be the whole-space weighted induced-field stiffness, with exact stationary
vacuum rows. On material nodes define the positive consistent weighted mass
matrix `M_ij=integral_plasma r³ N_i N_j dr dz`. The weak current representation
is `Jtheta=r q`, with `M q=(K u)_plasma/mu0`, u=a_ind. This includes interface
normal-derivative mismatch in the finite-element weak residual; it is not an
unresolved physical delta-function surface-current closure.

For constant prescribed resistivity eta>=0 and fixed uniform bias,

```
M a_material_dot = -(eta/mu0) (K u)_plasma,
Qdot = 2 pi eta q^T M q >= 0,
dE_ind/dt = -Qdot.
```

An implicit midpoint solve pairs magnetic loss with positive cell heating
`Q_c=2 pi eta dt integral_cell r³ q_mid² dr dz`. Add this physical heat to
cell internal energy and reconstruct the cell entropy invariant at fixed
volume. Do not double count it or change kinetic energy to repair a residual.
With a changing uniform bias, `u_dot=a_material_dot-Bdot/2`; count prescribed
source work using the same midpoint magnetic load and the actual bias increment.
Vacuum nodes still have zero material mass/resistivity and satisfy elliptic rows.

For varying eta, define `R_ij=integral_plasma eta r³ N_i N_j dr dz`, then
`M a_material_dot=-R q` and `Qdot=2pi q^T R q`. Do not replace a noncommuting
consistent mass/resistivity pairing with an unproved nodal multiplication.
Representation of current, interface regularity and local heat allocation must
be compared to continuum profiles independently of the global energy identity.

For constant eta the coupled material/vacuum midpoint matrix can use M on
material rows and zeros elsewhere, plus `dt eta K/(2mu0)`. This avoids explicitly
forming a dense vacuum Schur complement. The eta=0 limit needs separate handling
because zero vacuum rows in M alone cannot determine the exterior state.

## Independent spherical magnetic-diffusion eigenmode

In a stationary conducting sphere of radius R, use a=a(s,t), s²=r²+z², with
no imposed field, constant magnetic diffusivity D=eta/mu0 and true vacuum outside.
The induction equation and matching conditions are

```
a_t = D (a_ss + 4 a_s/s),
a_s(0)=0,  a_s(R)=-3 a(R)/R,
a_outside(s)=a(R) (R/s)^3.
```

For k=n pi/R and central initial a=A, a regular exact solution is

```
a_inside(s,t) = 3 A j1(k s)/(k s) exp(-D k² t),
a(0,t) = A exp(-D k² t),
Jtheta_inside = (k²/mu0) r a_inside; Jtheta_outside=0.
```

The Robin condition follows because
`d/ds[j1(ks)/(ks)] + 3 j1(ks)/(k s²) = sin(ks)/(k s²)`;
it vanishes at kR=n pi. Field and normal derivative match the vacuum dipole,
so no surface-current sheet is part of this continuum solution.

Using `integral_0^R s^4 a(s,0)^2 ds=9 A² R/(2k^4)`, independent whole-space
energy and Ohmic-power targets are

```
E_ind(t) = 6 pi A² R/(mu0 k²) exp(-2D k² t),
P_ohmic(t) = 2D k² E_ind(t),
Q(0,t) = E_ind(0) [1-exp(-2D k² t)].
```

For R=1 m, A=sqrt(mu0) T, n=1 and D=0.1 m²/s, E(0)=6/pi J and
P(0)=1.2 pi W. A separate 100-point radial Gauss integration of the gradient
energy plus analytic infinite exterior gave 1.909859317102750 J versus
1.909859317102744 J analytically; direct J² integration gave 3.769911184307752 W.
This checks the derived continuum formulas, not a numerical induction solver.

Required numerical studies: radial/angular and timestep refinement, current and
field profiles including the boundary, exact decay rate, local and integrated
heat, matching/mode/quadrature controls, regular-axis limits and atomic rejected
states. The frozen-geometry fixture is not a freely evolving force balance:
its geometry must be explicitly held fixed for this induction component test.
Then verify coupled ideal/resistive time splitting, moving-source work and
mesh sensitivity before interpreting any dynamic experiment. Material transport
and Hall applicability remain separate physical requirements.
