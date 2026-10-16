# Atomic reactions on actual material-support overlaps

`MaterialReactionOverlap(source, reference_density, quadrature=3)` gathers the
actual `MaterialSpecies` state onto Gauss points on the union of its physical
cell edges and remote field ends. Each material fluid has an explicit atomic
level/electron index; masses and charges must match exactly. There is no common
support, quasineutral projection, floor, or assigned population for empty levels.

`from_material` preserves the represented number, linear velocity and quadratic
kinetic energy. Two or more Gauss points integrate these initial extensive
moments exactly. The inherited source uses one fixed reference density, without
renormalizing cells to unit nucleus number. Empty exterior rows remain exactly
zero. Existing fixed-frame coordinates are retained; laboratory velocities and
moments are restored in outputs without discarding actual charge residuals.
MaterialSpeciesSet is also accepted: atomic IDs can repeat for disconnected
blocks, but two occupied blocks of the same fluid at one point reject. An empty
block never overwrites an occupied block's moments during gathering.

`step` advances the existing nonlinear exponential-midpoint atomic moment source
at every physical quadrature point. New levels/ions/electrons may appear where
reactants overlap, even beyond the former support of a product species. It keeps
the full local product moments, binding and photon reservoirs together. No
interpolation of separate charge products onto old meshes is performed. This
therefore retains pointwise reaction charge conservation, which is necessary
for the unchanged Gauss field during a reaction-only substep. Finite quadrature
of the nonlinear reaction solution is an explicit spatial approximation.

`fields` exposes physical number, velocity, temperature, charge and current at
the stored points. `totals` independently integrates number, mass, vector
momentum, translational + binding + produced photon energy, plus the Gauss
field through the exactly empty exterior. This budget is observed, never imposed.
The field diagnostic integrates interval-average charge; charge was initially
constant on each overlap and remains so under local reactions up to roundoff.
No physical force, field work, photon transport or motion occurs in this source
substep. Existing small-drift Maxwellian, capture, missing-channel and photon
recoil assumptions remain those of atomic_moments.md and the supplied network.

**Integration boundary:** the result is `ReactionOverlapState`, not an advanced
`SpeciesMaterialState`. It deliberately retains unprojected local moments.
Returning them to a moving representation must preserve local birth charge,
variable particle amounts, momentum/energy and newly occupied support; this
remaining transport/reaction integration is not implemented by freezing numbers,
discarding products, distributing heat to close a budget or rescaling old J.
The full moving reactive goal remains open; this is not formation adoption.
`reactive_material` now supplies a first explicit moving/source composition.
Its cold translating neutral decay limit completes, while the full hot charged
argon attempt exposes mesh/timestep collapse; see ATOMIC_MOVING_FIELDS.md.
