# Analytic spherical FRC reference

Independent continuum reference for a future moving plasma/vacuum boundary.
It prescribes an ideal, static poloidal Grad–Shafranov equilibrium, with a
spherical separatrix and a uniform-bias-plus-dipole vacuum extending to infinity.
It does not advance the existing fixed-chamber solver or claim formation.

`fields(r,z)` returns SI psi, Br, Bz, pressure and toroidal current, including
regular axis limits. Interior and exterior magnetic fields are continuous at
the sphere; normal B and pressure vanish there, so total traction matches
without a sheet current or exterior gas. Volume current jumps at the boundary.
`integrals(gamma)` gives independently derived continuum flux/current/energy
targets. Uniform-bias self energy is infinite and explicitly excluded from
the finite induced-field energy. No density or temperature is silently chosen.

See `docs/FREE_BOUNDARY_REFERENCE.md` for equations, sources and limitations.
Independent derivative, interface and quadrature tests verify this analytic
reference. It remains separate from evolving numerical solvers and their
current fixed-chamber source fingerprints.
