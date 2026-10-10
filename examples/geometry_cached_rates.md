# Reuse original field geometry within a complete batch

CompatibleGeometry contains state-independent field bases, physical quadrature,
mass matrices/factors, magnetic Hessian and original Faraday incidence. Within
one run, retain it for exactly matching cylinder/field mesh, polynomial degree
and quadrature settings. All supplied option values and types enter the key.
No solver matrix is approximated or re-factored by another algorithm.

Each case still constructs its actual fluid space, field mapping, state,
density, wall geometry, current inertia, species/transport and induction
operators. StrongWallGeometry copies the original geometry attributes and
constructs its own wall projection/curl; it does not mutate the original.
Only the already computed immutable original geometry is shared. No density,
initial state, forcing, solution or rate is shared between cases.

The cache starts empty for every full run and full replay. Summary metadata
records constructions and hits. Full source closure includes this serial
dispatch and the original constructor; bindings restore on all exits.
Run parallel variants in separate processes.

Verify with complete original/new signed workflows on distinct fluid grids
sharing one field grid, requiring every input/rate, solver diagnostic, physical
fingerprint and norm to agree bitwise; then fully replay the new batch. This
is an end-to-end workflow, with no unit tests. Cached geometry only changes
setup reuse, not equations, boundaries, residuals, diffusion or physics.
