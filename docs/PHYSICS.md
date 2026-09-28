# Physics contract

All internal quantities are SI. Cylindrical coordinates are right handed
`(r, theta, z)`; positive loop current flows in +theta and produces +Bz on axis.
The axisymmetric domain includes `r=0`, with no theta dependence.

## Progressive models

| Stage | Status / equations | Missing physics |
|---|---|---|
| 0 | Implemented magnetoquasistatic vacuum coil superposition | Radiation, displacement current, conducting walls |
| 1 | Implemented stationary scalar-conductivity induction and surrounding vacuum | Motion, thermodynamic feedback, transport closure |
| 2 | Conservative bounded MHD reference plus verified smooth fixed-chamber gas/coil/vacuum coupling | Free plasma interfaces, physical ends, magnetic shocks and formation workflow remain incomplete |
| 3 | Research: energy, transport, ionization, neutrals, radiation | Atomic data and closure choices unresolved |
| 4 | Research: Hall / two-fluid | Electron pressure and inertia closures unresolved |
| 5 | Research: 3-D / hybrid / kinetic | Nonaxisymmetric modes and particle confinement |

## Flux convention and Stage 0

`psi = r A_theta = axial flux / (2 pi)` [Wb, conventionally Wb/rad].
`Br = -psi_z/r`, `Bz = psi_r/r`, `Etheta = -psi_t/r`.
The axis gauge is `psi(0,z)=0`. Smooth fields have
`psi = Bz(0,z) r^2/2 + O(r^4)`, `Br(0,z)=Etheta(0,z)=0`.
Uniform bias contributes `psi=Bbias*r^2/2`. It is a prescribed static background.

For a loop radius a, axial offset h=z-zc and ampere-turn current NI, define
`D=(a+r)^2+h^2`, `q=(a-r)^2+h^2`, `m=4ar/D`.
With complete elliptic integrals using the **parameter m**, not modulus k:

```
psi = mu0 NI sqrt(D)/(2 pi) [(1-m/2) K(m) - E(m)]
Br  = mu0 NI h/(2 pi r sqrt(D)) [-K(m)+(a^2+r^2+h^2) E(m)/q]
Bz  = mu0 NI /(2 pi sqrt(D)) [K(m)+(a^2-r^2-h^2) E(m)/q]
```

The axis limit is `Bz=mu0 NI a^2/[2(a^2+h^2)^(3/2)]`.
Small-m series for psi and near-axis field series avoid cancellation.
Filament fields are singular on the wire: evaluating sufficiently close to a
filament raises an error, rather than clipping or inventing a wire radius.
The analytic field extends to infinity; no mesh boundary is imposed in Stage 0.
Currents are prescribed functions; their derivatives determine Etheta.
Sampled currents use linear interpolation and piecewise slopes (right derivative
at interior knots, zero beyond the table); callback derivatives must be supplied.
Finite pulses use smooth endpoint ramps. Radiation is neglected: require
`L/(c*t_wave) << 1`. A waveform is a source specification, not a circuit design.

## Stage 1: stationary conductor coupled to vacuum

Let `u=psi_induced`, `psi=psi_external+u`, with prescribed time-independent
conductivity sigma(r,z) inside a cylindrical plasma. Outside it sigma is exactly
zero, including the quartz wall approximation (mu=mu0; dielectric response absent).
There is no artificial conducting background and no wall surface current.

```
Delta* = d_rr - (1/r)d_r + d_zz
Jtheta = -Delta* u/(mu0 r)
sigma (u_t + psi_external,t) = Delta* u / mu0    in conductor
Delta* u = 0                                   in vacuum
```

These follow from Faraday, Ampere without displacement current, and `J=sigma E`
at v=0. The externally prescribed coil field is analytically source-free in
the plasma. It is not passed through a discrete Laplacian to manufacture J.
Induced flux and its normal derivative are coupled through one elliptic mesh;
sigma is node sampled, so curved/interface positions are grid approximations.

Boundary conditions: u=0 on axis (regularity), radial outer boundary and both
axial outer boundaries. The latter approximate decaying induced fields at a
finite distance, **not a physical quartz-wall boundary**. Domain enlargement is
required. Coils may lie in the vacuum mesh; only their analytic external field
is used. Initial u=0 means a plasma initially penetrated by the initial field,
with no pre-existing induced current. Arbitrary initial u must solve the vacuum
constraint; the solver checks it.

Discretization: uniform nodal (r,z) mesh, interior unknowns. A symmetric flux
stiffness K discretizes `-div(grad(u)/r)` with edge weights
`dz/(r_face dr)` and `dr/(r_node dz)`. The diagonal mass is
`M=mu0 sigma dr dz/r`. Thus `M udot+K u=-M psi_external,dot`.
Backward Euler uses exact external-flux increments:
`(M+dt K) u_new = M u_old - M (psi_ext,new-psi_ext,old)`.
Vacuum rows have zero mass and remain elliptic. Sparse LU is cached by dt.
This is first order in time and nominally second order in smooth space; interface
stair steps can reduce observed order. No smoothing, floors, or hyperdiffusion.

Etheta and Jtheta at each step are backward-interval quantities; current obeys
Ohm's law and the discrete Ampere equation. Reconstructed nodal B is for plots
and probes; face-integrated fluxes from psi telescope to zero cell divergence
to roundoff. Finite-difference divergence of interpolated nodal B need not be zero.

## Energy accounting

The exact discrete induced-field energy is `Uind=pi/mu0 * u.T K u` (includes
vacuum within the outer domain). Dissipated heat per step is
`Q=dt sum(Jtheta^2/sigma * 2pi r dr dz)` over conductor nodes.
External work on induced currents is `W=-sum(2pi Jtheta Delta psi_ext dr dz)`.
Backward Euler removes `Dnum=pi/mu0 Delta u.T K Delta u >=0`.
The verified identity is `Delta Uind+Q+Dnum-W=0` to solver precision.
Q is accumulated but does not update temperature or sigma. Magnetic energy from
total B is separately integrated over the chamber; this is not a closed-system
energy budget and excludes source circuits and energy outside the chamber.

## Diagnostics and interpretation

Field reversal is measured relative to the sign of nonzero Bbias within the
plasma; on-axis reversal is separate. O/X points are numerical critical-point
candidates of the interpolated flux, classified by the Hessian. Axis nulls are
reported separately. Closed contours must close inside the domain and enclose
a detected O point. Sampled contours are evidence at the mesh resolution, not a
proof that every island is detected. Separatrix estimates require suitable axis
nulls or X points and are flagged as estimates; ambiguous cases return null.
Lifetime is a sampled lower bound; brief events between outputs can be missed.
Excluded-flux radius is not inferred from reversal; its experimental definition
remains open. Stage 1 cannot establish force balance, an MHD equilibrium,
stability, a successful FRC, or particle confinement.

Characteristic scales use explicitly supplied species mass, charge, ion density,
electron/ion temperatures and a stated reference |Bbias| and plasma radius L:
`p=e(ni Ti+ne Te)`, `beta=2mu0 p/B^2`, `vA=B/sqrt(mu0 ni mi)`,
`tauA=L/vA`, `taud=mu0 sigma L^2`, `S=mu0 sigma L vA`,
`Rm=mu0 sigma L U` (U=0 in Stage 1), `di=sqrt(mi/(mu0 ni Z^2 e^2))`,
`omega_cs=|qs|B/ms`, `rho_s=sqrt(2e Ts/ms)/omega_cs`.
These are reference estimates, not local guarantees near magnetic nulls.
An effective Drude collision frequency `ne e^2/(me sigma)` is labeled an
inference from prescribed sigma, not a validated plasma collision model.
Physical collision frequencies, mean free paths and neutral effects remain
unknown until a transport model is chosen. Hall and gyroradius warnings are
heuristics for missing physics, not automatic permission to apply a fluid model.

## Sources

- Simpson et al., NASA, circular-loop field expressions:
  https://ntrs.nasa.gov/citations/20010038494
- SciPy elliptic-integral parameter convention:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipk.html
- Fitzpatrick, induction equation and diffusion:
  https://farside.ph.utexas.edu/teaching/plasma/lectures1/node123.html
- NRL Plasma Formulary, characteristic scales (converted to SI here):
  https://www.nrl.navy.mil/Portals/38/PDF%20Files/NRL_Formulary_2019.pdf

Numerical discretization and its discrete energy identity are derived above;
verification measurements are generated in docs/VERIFICATION.md.

The new dynamic kernel's equations, boundary restrictions and explicit numerical
diffusion are in [MHD_NUMERICS.md](MHD_NUMERICS.md); measured wave, shock,
advection, resistive-heating and time-convergence checks are in
[MHD_VERIFICATION.md](MHD_VERIFICATION.md). Its gas-filled conducting-wall
verification cylinder must not be confused with the quartz/plasma/vacuum setup.
The [independent closed-flux equilibrium checks](EQUILIBRIUM_VERIFICATION.md)
test curved magnetic stress, reversed Bz and residual artificial motion. They
start from a prescribed equilibrium and do not demonstrate dynamic formation.

The separate [dynamic vacuum coupling](COUPLED_NUMERICS.md) uses gas-only energy
and global induced self-field energy with an adjoint force/induction work pair.
It includes moving gas, scalar resistivity and prescribed coil work without an
exterior fluid atmosphere. [Independent measurements](COUPLED_VERIFICATION.md)
show smooth convergence and quantify the finite magnetic boundary. Its fixed,
gas-filled reflecting chamber is still a verification restriction, not a moving
free plasma or validated quartz/end model. The response-energy budget excludes
external-field self/mutual energy; it is not total electrical input from a coil.
The optional [infinite magnetic exterior](OPEN_BOUNDARY.md) now replaces the
fixed far-field condition with decaying vacuum multipoles. Its field energy
beyond the grid is included in the same response budget, and [independent
verification](OPEN_BOUNDARY_VERIFICATION.md) separates mode truncation, spatial
error and matching-surface sensitivity. This changes the magnetic boundary only;
it does not create an evolving free fluid surface or establish transport validity.
