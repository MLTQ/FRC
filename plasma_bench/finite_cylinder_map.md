# Finite fluid/field physical maps

FiniteCylinderMap retains the original actual regular radial field definitions,
force-work adjoint, velocity quotient and density pullback. Its axial map is the
finite Lobatto nodal polynomial, evaluated through a well-conditioned Legendre
basis. It is not a periodic interpolant. FiniteFluidDifferential differentiates
that same physical polynomial and reuses the cylindrical gradient/stress pair.

FiniteSpecies retains the original skew two-species entropy/current transport
and compatible potential split. The potential shape and gradient adjoint now
include both distinct axial endpoint nodes. No species or equation is removed.
All actual point/quadrature densities remain positive with no floors. This map
requires FiniteCylinderFluid and FiniteCylinderGeometry; mixed periodic/finite
representations are not supported.

Seven full24us coil pulses now exercise this coupled path. Independent limits,
spatial/time comparisons, transport flags and remaining exterior-resolution
limits are recorded in docs/FINITE_CYLINDER_RESULTS.md.
