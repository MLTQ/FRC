# Actual moving species, pressure and charge separation

Runs three physical equation sectors on 32/64/128 moving cells, with every
accepted geometry and integrated species state preserved. Reactions/collisions
are explicitly disabled in this initial moving-operator study; species masses,
momenta and total energies are physical. The inert neutral/one-ion layout is not
claimed to be a complete argon atomic model. No unit-test suite is used.

`--mesh-motion affine` (default) distributes coordinate motion between the
tracked fronts so outer cells do not grow while the interior stays fixed.
Physical species velocities remain independent, not constrained to affine flow.
`--mesh-motion barycentric` retains the original material-weighted coordinate
choice. Both solve the same ALE conservation laws and retain all exterior energy.
Cell-width ratios are reported to expose distortion, not to repair a mesh.
The current conservative reconstruction preserves each cell's mean total energy
as well as its number/momentum; admissibility limiting affects only face slopes.
The reported limited fraction counts all occupied-species reconstruction calls,
including timestep queries. Earlier primitive-reconstruction failures remain
archived and are not overwritten.

## Independent neutral reference

Neutral argon starts at rest with half-width a=0.025 m, T=0.1 eV and Nref=1e19/m³.
The duration a/(4c0) precedes interaction of the two rarefaction heads. The exact
Euler vacuum-fan solution is given in atomic_expansion.md. Reference integrals
split at heads and fronts and use polynomial-exact Gaussian quadrature. Errors
include any analytic matter outside the numerical front. Exact emptiness outside
the tracked front is reported separately from the front's position error.

## Independent charged warm affine expansion

A positive-ion slab initially has n=n0=1e12/m³ for |z|<a=0.025 m, no electrons,
and T=T0[1-(z/a)²], T0=1 eV. The remote electric boundaries are at ±0.2 m, with
E_left=-e n0 a/epsilon0 and the right field determined from actual total charge.
There is no fictitious neutralizing volume charge. Fixed remote surfaces supply
the electric boundary stress, whose net force is zero by symmetry.

With gamma=5/3, the continuum Euler–Poisson solution remains affine:

    z=s(t) X, n=n0/s, v=s_dot X,
    T=T0 s^(1-gamma)[1-(X/a)²],
    E=e n0 X/epsilon0 inside; constant at either exterior side,
    s_ddot=A s^-gamma + omega_pi²,
    A=2 e T0/(m_i a²), omega_pi²=e² n0/(epsilon0 m_i).

The duration to s=1.5 is obtained from the independent first integral
s_dot²=2A[1-s^(1-gamma)]/(gamma-1)+2 omega_pi²(s-1), by scalar quadrature with
s=1+y² to remove the initial square-root singularity. This does not integrate
the numerical moving-volume equations. Initial and final cell integrals retain
the exact quadratic pressure and kinetic profiles. Exterior field energy is
included in the physical conservation budget.

## Neutralized electron–ion free expansion

A 0.02 m-wide slab at n_i=n_e=1e14/m³ and Ti/Te=0.2/2 eV evolves for 20 ns.
All electric forces, relative species fluxes, charge separation and boundary
motion evolve; neutrality is only an initial condition. The outermost Euler
vacuum fan tracks common support, while each species may spread differently
inside it. This is a nonlinear fluid exploration, not a kinetic sheath or
laboratory-formation reference. Final field/density/front differences are measured
between grids, without treating one grid as exact truth.

The operator's acoustic, mesh-crossing and plasma-frequency restrictions select
actual steps. No failed state is clipped, repaired or retried at hidden settings.
Failures retain every accepted prefix and attempted-step diagnosis. Scientific
figures and all source identities are saved. Transverse fields, cylindrical
geometry, actual reactions and physical wall/open-end transport remain separate.

Current results are in `docs/ATOMIC_MOVING_FIELDS.md`. The finest plasma cases
stop before the requested 20 ns and cannot be compared as equal-time endpoints.
Figures now label the actual endpoint time and completion status; the first
archived plot is preserved and separately annotated from its saved arrays.
