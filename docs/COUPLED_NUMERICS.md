# Developmental dynamic plasma–vacuum coupling

## Scope and boundaries

`VacuumCoupledMHD` is a second dynamic reference, alongside the conservative
bounded `AxisymmetricMHD`. Fluid exists only in a declared **gas-filled, fixed
cylindrical chamber**. Radial walls and axial ends reflect normal gas velocity
and pass magnetic flux. The larger exterior has no density, pressure or
conductivity: only a quasistatic elliptic magnetic field. The finite outer
boundary sets induced psi to zero. Prescribed analytic filament coils and a
uniform bias enter through the independent `VacuumField` interface.
Prescribed source filaments must lie outside the fluid chamber. The exterior
zero-current constraint concerns the induced response, excluding those sources.

This is not yet a detached free plasma surface, an open-ended tube, dielectric
wall physics, or physical transport for a selected gas. No numerical atmosphere
or artificial exterior conductivity is introduced. The main experimental
`stage=2` gate remains in place. Displacement current, toroidal magnetic field,
rotation, Hall terms, ionization and kinetic effects are absent.

## Equations and energy convention

Use SI and psi=magnetic flux/(2 pi), Br=-psi_z/r, Bz=psi_r/r. Write total flux as
psi=u+psi_external. The existing symmetric stiffness K discretizes
-div(grad(u)/r); its nodal integral Ku approximates mu0 Jtheta dA. The induced
self-field energy is W=pi/mu0 u^T K u over the entire magnetic domain. Gas energy
is internal plus kinetic; it does **not** contain a second magnetic-energy term.

Let P mark plasma nodes, V vacuum nodes, A_i be their dual patch area inside the
actual gas cylinder, M_v the diagonal matrix of cell masses repeated for both
velocity components, and R_i=eta r_i/(mu0 A_i) on P and zero on V. At chamber
edges/corners, A_i is the appropriate half/quarter patch. Define C(psi) v as the
nodal interpolation of v dot grad(psi), with zero normal velocity at each wall.
For the magnetic substep density is fixed and

```
P [u_dot + C(psi) v + psi_external_dot + R K u] = 0,
(K u)_V = 0,
M_v v_dot = C(psi)^T g,          g = (2 pi/mu0) K u.
```

The paired force is a weak discretization of Jtheta cross B_total: in the
continuum Jtheta v dot grad(psi)/r equals v dot (J cross B). The transpose is
essential for discrete work cancellation; it does not by itself ensure accurate
pointwise force, local momentum fluxes, shock jumps or equilibrium preservation.
Those need independent verification. Centroid interpolation is linear and may
extrapolate tangential velocity at boundaries; it is not a positivity limiter.

## Time advancement and nonlinear solve

Strang splitting applies a gas-only Euler half-step, a magnetic midpoint full-step,
and an Euler half-step. Euler uses annular finite volumes, limited second-order
reconstruction, Rusanov fluxes and SSP RK2 with reflected gas ghosts. Its numerical
diffusion acts on gas transport. The magnetic step uses centered derivatives and
implicit midpoint, with no added magnetic artificial resistivity. This smooth
energy-stable formulation has no established magnetic-shock monotonicity guarantee.

At tau=dt/2, eliminating midpoint velocity gives, for a fixed iterate of C,

```
H = R + tau (2 pi/mu0) C M_v^-1 C^T,
[P + tau H K + V K] u_mid
    = P (u_old - delta(psi_external)/2) - tau C v_old.
v_mid = v_old + tau M_v^-1 C^T g_mid.
```

Here P,V in matrix expressions denote row-selector diagonal matrices. C has
zero vacuum rows. Sparse direct solves and fixed-point iteration update C from
u_mid plus the average external flux. Both iterate change and the nonlinear
induction residual must meet the tolerance. Nonconvergence raises; no partial
state, clock or cumulative energy is committed. No cell density/pressure floor
or failed-state energy reset is used. A conservative fast-wave CFL and explicit
diffusion accuracy bound still cap dt, although the magnetic step is implicit.
An optional `max_timestep` resolves known waveform features. A fluid CFL cannot
discover an arbitrary narrow callback pulse; source time resolution must be
declared and tested. The provided coil fixture caps dt at one twentieth of rise.

## Measured work identity

Midpoint yields delta W = g_mid^T delta u. The velocity update gives exact
discrete kinetic work v_mid^T M_v delta v. Positive Ohmic heat is

```
Q = dt (2 pi/mu0) sum_i R_i (K u_mid)_i^2 >= 0,
W_source = -g_mid,P^T delta(psi_external)_P,
delta(W + E_gas) = W_source.
```

Heat is distributed to adjoining fluid cells with nonnegative fractions summing
to one. Gas energy changes by the actual kinetic work plus Q in this physical
subsystem; internal energy is constant when eta=0. The separate gas transport
step conserves mass and gas energy under closed walls. Diagnostics report the
global budget, magnetic exchange residual, vacuum current residual, nonlinear
residual and compatible face divergence.

This **response-energy** identity includes induced self-field energy and gas
energy. It excludes external-field self energy and mutual energy. W_source is
the applied electric field's work on plasma current, not an estimate of complete
coil-circuit input energy or a basis for equal-delivered-energy hardware claims.

## Relation to published methods and verification requirements

Splitting compressible gas dynamics from magnetic work is supported by
[Dao, Nazarov and Tomas](https://arxiv.org/abs/2310.18467); their finite-element
method includes invariant-domain and energy analysis. The Cartesian
[Pang and Wu PPCT method](https://arxiv.org/abs/2410.05173) uses an implicit
energy-conserving magnetic subsystem with Strang splitting. The present annular,
fixed-interface vacuum coupling is independently derived. It does not inherit
either paper's full positivity or convergence theorems.

Required checks include continuum cell-averaged J cross B, independent translating
flux and Gaussian diffusion solutions, temporal refinement, coil work and Joule
heat, and an infinite-space circular-loop Green-function comparison of the
finite magnetic domain. Passing these smooth checks does not establish physical
formation accuracy. Free interfaces, end effects, magnetic shocks, strong
compression, transport validity and converged formation remain open requirements.
