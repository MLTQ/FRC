# Axisymmetric dynamic MHD: governing and numerical contract

This developmental kernel evolves a compressible conducting fluid throughout a
cylinder. It is a **verification environment**, not yet the coil-driven plasma
surrounded by vacuum. The experimental Stage 2 entry point remains gated until
that coupling is implemented and tested. No vacuum density or conductivity
floor is introduced to obtain an attractive formation plot.

## Equations and units

SI; ideal-gas gamma, scalar constant resistivity eta [ohm m]. No theta dependence,
no azimuthal velocity and no toroidal magnetic field. These restrictions retain
poloidal FRC dynamics but omit rotation, Hall-generated components and all
nonaxisymmetric instability modes. Energy density is
`E=p/(gamma-1)+rho*(vr^2+vz^2)/2+(Br^2+Bz^2)/(2 mu0)`.

Define `pt=p+B^2/(2 mu0)`, `U=(rho,rho vr,rho vz,E)` and

```
Fr = (rho vr, rho vr^2+pt-Br^2/mu0,
      rho vr vz-Br Bz/mu0, (E+pt)vr-(v.B)Br/mu0)
Fz = (rho vz, rho vr vz-Br Bz/mu0,
      rho vz^2+pt-Bz^2/mu0, (E+pt)vz-(v.B)Bz/mu0)
dU/dt = -[d(r Fr)/dr]/r - dFz/dz + (0,pt/r,0,0)
psi_t = -r Etheta
Etheta = vr Bz-vz Br+eta Jtheta
Jtheta = -Delta* psi/(mu0 r)
```

Resistive energy flux adds `(eta Jtheta Bz/mu0)` radially and
`(-eta Jtheta Br/mu0)` axially. There is no separate Ohmic source in **total**
energy; magnetic energy becomes internal energy. Adding eta J² to total energy
as well would double count heating. Physical heat is diagnosed separately.

## Mesh and boundaries

Cell-centered volume-averaged conserved variables have shape `(nr,nz,4)`.
Flux psi lives at `(nr+1,nz+1)` cell corners. Radial faces include r=0;
axial endpoints are periodic duplicates. Cell volumes are exact annular
volumes `pi*(r_outer^2-r_inner^2)*dz`. Magnetic face averages are

```
Br(i,j) = -(psi(i,j+1)-psi(i,j))/(r_face(i)*dz)
Bz(i,j) = 2*(psi(i+1,j)-psi(i,j))/(r_outer^2-r_inner^2)
```

The axis has zero radial area, Br=0 and psi=0. The cell face flux divergence
telescopes exactly. The outer cylindrical boundary is impermeable, perfectly
conducting with constant psi along z (Br=0) and tangential fluid slip. The outer
psi value stays fixed. No axial end loss exists with periodic z. These are
declared verification boundaries, not quartz/vacuum boundary conditions.

## Discretization

The initial reference method uses first-order piecewise-constant states and
local Lax–Friedrichs (Rusanov) face fluxes. Both face states use the single
staggered normal magnetic field, preserving their thermal pressure when their
face total energy is formed. This supplies conservative finite-volume mass,
momentum and total-energy updates. Exact annular face areas and the geometric
radial source balance a uniform stationary pressure/axial field at the axis.

The ideal corner EMF uses four neighboring fluid states, staggered transverse
field states and **directional upwind magnetic jumps**. In axisymmetric signs:
`Etheta = mean(vr Bz-vz Br) - ar/2*(Bz_R-Bz_L)
                              + az/2*(Br_U-Br_D)`.
Here ar, az bound directional fast-wave speeds in neighboring states. This
retains the full Rusanov dissipation in a one-dimensional limit. A simple
average of four face electric fields would halve it, so is not used.

The resistive corner EMF uses second-order centered Delta* and the regular axis
limit. Its face-averaged Poynting flux enters total energy. The method is explicit;
its time step must satisfy both fast-wave and resistive diffusion restrictions.
SSP RK2 advances fluid and flux together. The spatial hyperbolic scheme remains
first order; temporal order is two for the semidiscrete system.

## Stabilization and conservation

The Rusanov/UCT terms introduce numerical viscosity and magnetic diffusion of
order wave-speed times grid spacing. They are not physical resistivity, cannot
be adjusted to force FRC formation, and must be reduced/measured by refinement.
No density/pressure floors, magnetic-energy resets or silent positivity repairs
are allowed. Invalid states raise, identifying the failed quantity; an excessive
requested step is rejected. Diagnostic budgets separately retain total,
thermal, kinetic and magnetic energies and physical eta J² estimates.

Global mass and total energy are conserved to roundoff under the stated closed
radial and periodic axial boundaries. Conservation alone does not prove correct
energy partition, wave propagation, topology or transport. These need independent
analytic/manufactured benchmarks and measured refinement. A closed-flux initial
condition, if used as a verification fixture, is not spontaneous FRC formation.

## Method references

- [Stone et al. (2008), Athena](https://arxiv.org/abs/0804.0402): conservative
  fluid quantities with staggered constrained magnetic transport.
- [Skinner & Ostriker (2010), cylindrical MHD](https://arxiv.org/abs/1004.2487):
  cylindrical finite-volume geometry and source terms.
- [Mignone & Del Zanna (2020), UCT](https://arxiv.org/abs/2004.10542): consistent
  upwind corner EMFs and correct one-dimensional reduction.

The implementation here is an independently written restricted reference kernel,
not Athena or a claimed reproduction of its higher-order solver.
