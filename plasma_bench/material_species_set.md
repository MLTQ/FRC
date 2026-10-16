# Different material meshes for different occupied species

`MaterialSpeciesSet(models, species_indices)` advances occupied one-species
MaterialSpecies blocks, each with its own count, physical faces, particle amounts,
three nodal velocities and cell thermal energies. `MaterialSetState` stores their
states without padding vacuum with cells. Atomic IDs permit disconnected support
pieces; later overlap of pieces belonging to one species rejects rather than
silently adopting a multistream closure. Every block shares the remote field ends.

The pressure, numerical stress, consistent mass and work are exactly those in
species_material.md. The common electrostatic field is integrated from all actual
charge on the union of their edges. Each block's existing two-point electric
force integral uses this common field, not its isolated self-field. RK4 advances
all blocks together. Acoustic restrictions and the sum-of-density plasma-rate
bound retain the physical electron mass. No floor, remap or energy repair.

`sample` returns one row per support block and the total J/E. `totals` computes
actual matter energy/momentum and the shared interior/empty-exterior field energy
once, avoiding double-counting. `faces` exposes all physical faces for the
reaction-overlap gather; unlike rectangular MaterialSpecies, block lengths vary.

An explicit fixed `reference_velocity` stores nodal velocities relative to that
frame; geometry advances with v_relative+V_frame. Forces and thermal work use
relative velocity differences. Sampling and total energy/momentum restore lab
motion, including charge-current transport. This coordinate choice avoids
subtracting large common bulk energy from an exactly cold translating gas;
it does not reset temperature, velocity, or charge.

Numbers remain fixed within each transport substep. A reaction step must supply
new physical blocks/amounts from its complete moments; this class alone does not
implement reaction deposition, mesh coarsening, collisions, radial geometry or
physical walls. Source-to-transport projection is a separate numerical operator
whose local charge, moment/energy error and admissibility must be assessed.

Opt-in MaterialSubcells blocks have fine density geometry and fewer velocity
nodes. The set delegates their compatible force/work and acoustic restriction
and uses their represented fine endpoint velocities in kinetic/momentum totals.
Every density cut still participates in the shared Gauss field/plasma bound.
Scalar E sampling converts extended geometry to float64 for np.interp output;
force integration retains the original extended coordinates and field arrays.
Subcell thermal coordinates are divided by a fixed per-block particle scale;
the set restores their physical value for the energy budget.
Subcell number density and particle totals are likewise recovered from factored
number coordinates; fine physical amounts can be below float64 SI range without
being deleted from the internal state. Totals multiply the scale after summing.
