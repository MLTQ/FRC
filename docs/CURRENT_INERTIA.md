# Current inertia before a complete plasma-edge model

The sharp heated-edge audit in SURFACE_COMPATIBILITY remains unresolved. This
component adds an independently testable term needed when carrier response is
not instantaneous. It does not claim that this term alone repairs that surface.

## Equations and restricted physical limit

At fixed density and fixed geometry, the toroidal sector is

```
curl B = mu0 J,       div B = 0,       d_t B = -curl E,
Etheta = eta Jtheta + [me/(ne e²)] d_t Jtheta,
de² = me/(mu0 ne e²),   eta >= 0.
```

Electron inertia and current energy follow from electron momentum. The full
generalized Ohm law also has Hall, pressure and convective terms: see Martin,
[equations 2.11–2.18 and section 2.3](https://ecommons.cornell.edu/server/api/core/bitstreams/17448f20-41b4-4457-a472-3ff41ab87529/content).
That work studies different plasma configurations; it is not validation here.

This component isolates linear toroidal induction about rest and zero magnetic
field. Perturbation J and B are first order; J cross B, carrier self-advection
and fluid acceleration enter at higher order. A prescribed scalar relaxation
coefficient is part of this limit, not a complete frequency-dependent kinetic
conductivity. Density/temperature are fixed; positive collisional heat is
diagnosed rather than fed back. Finite bias, large drift or appreciable heating
requires additional terms and thermal evolution. The general moving solver
does not acquire this component until those equations and budgets are verified.

The exterior is exact magnetic vacuum with no carrier or material arrays.
The induced field decays at infinity; analytic uniform prescribed bias is
separate. Magnetoquasistatics still requires electromagnetic transit times small
relative to the evolution time. No displacement current or sheath is modeled.
An internal vacuum cell is rejected; it cannot be replaced by a density floor.

## Mixed finite-element equations

Use the existing a=psi/r² field, full-space stiffness K and plasma mass
`Mij=integral_plasma r³ Ni Nj dr dz`. Let `w=mu0 Jtheta/r`,
`Hij=integral r³ de² Ni Nj`, `Rij=integral r³ (eta/mu0) Ni Nj`.
Profiles are sampled at physical plasma quadrature points and then frozen.

```
K u = [M w; 0_vacuum],
M (d_t u_material + d_t a_external) + H d_t w + R w = 0.
```

The second equation is tested only in plasma. A sparse mixed midpoint solve
retains stationary vacuum rows and avoids dense Schur-complement assembly.
Solving for increments makes zero forcing give exactly zero update and reduces
cancellation when the physical change is small relative to the stored field.
Diagonal row/column equilibration reduces axis-conditioned roundoff in the
mixed solve; it is an algebraic scaling, not a physical regularization.
When both inertia and eta are exactly zero, the existing exact frozen-flux
branch preserves tiny material labels without subtracting and re-adding a
much larger bias. Its vacuum and Ampere residuals are still evaluated.
It allows explicitly prescribed positive density and nonnegative eta profiles.
The electron mass is its physical constant; zero inertia is a separately named
model comparison, not an artificial electron-mass tuning parameter.

Multiplying midpoint Ohm's law by `2 pi w_mid/mu0` gives

```
EB = (pi/mu0) u.T K u,
Ecurrent = (pi/mu0) w.T H w = integral me J²/(2 ne e²) dV,
Qcell = dt integral_cell eta J_mid² dV >= 0,
Wsource = -(2 pi/mu0) w_mid.T M Delta(a_external),
Delta(EB + Ecurrent) + sum(Qcell) = Wsource.
```

Energy is diagnosed independently, never repaired. Midpoint still permits
underresolved oscillatory decay, so time convergence is required. Variable
density changes current inertia, not resistivity; no eta taper suppresses heat.

## Independent continuum verification

For a sphere, `a=3 A j1(ks)/(ks) exp(-lambda t)`, `k=n pi/R`, matches
`a(R)(R/s)³` in vacuum. Substitution gives
`Jtheta=k² r a/mu0`, `lambda=D k²/(1+de² k²)`.
Its full-space magnetic energy is `6 pi A² R/(mu0 k²) exp(-2 lambda t)`;
current energy is `de² k²` times this, and Joule loss removes their sum.
Finite Maxwell/Faraday derivatives, direct volume integrals, mesh/time studies,
massless limits and nonuniform-coefficient work tests check distinct contracts.

These fixed-geometry benchmarks and conditional regime scans remain component
verification. They neither solve the heated material edge nor establish FRC
formation, ionization balance, particle confinement or three-dimensional stability.
