# Pressure-driven neutral expansion

Exercises nonzero pressure work, accelerating flow and a vacuum interface in the
new spatial species operator. All particles are neutral ground-state monatomic
argon with mass 39.948 u; reactions/collisions/fields are disabled explicitly.
The initial stationary slab has half-width a=0.025 m, T0=0.1 eV, Nref=1e19 m^-3
and exact zero matter outside, inside a 0.2 m periodic component domain.
The end time a/(4 c0) precedes interaction of the two rarefaction heads or physical
vacuum fronts with each other or the periodic boundary.

At the right edge, xi=(z-a)/t. The independent Euler solution is uniform for
xi<=-c0, exact vacuum for xi>=2c0/(gamma-1), and between them:

    u = 2(c0+xi)/(gamma+1)
    c/c0 = 2[1-(gamma-1)xi/(2c0)]/(gamma+1)
    n/n0 = (c/c0)^(2/(gamma-1))
    p/p0 = (c/c0)^(2 gamma/(gamma-1)).

The left edge is its reflection. Exact cell integrals split at both heads and
fronts before Gaussian quadrature; at gamma=5/3 the needed fan moments are
polynomials of degree at most five. The reference is independent of the numerical
flux. Grids 64/128/256 use fixed acoustic CFL, hence refine space and time together.
Every accepted state, analytic final state and number/momentum/energy budgets are
saved. Empty charged species remain exact zero.

Numerical front diffusion is measured by actual particle content and maximum
density in cells wholly beyond the exact vacuum front. It is not silently
replaced with an atmosphere, clipped away, or described as an exact free boundary.
This Eulerian operator is not adopted for the formation interface. The existing
moving plasma/vacuum boundary remains separate. This neutral Euler fixture is
not a collisionality claim for a laboratory gas or plasma edge.
