# Explicit source intervals and fine transport stepping

`SubcycledReactiveMaterial` uses source(dt/2), transport(dt), source(dt/2).
Transport(dt) is composed of actual existing RK4 material steps, each no larger
than transport_fraction times the current acoustic/plasma bound. That bound is
recomputed after reactions and after each accepted transport step. Physical
masses, all species, precise supports, fields, fine density cuts and actual U
remain unchanged in their meaning. There are no retries or forced large steps.

The caller explicitly chooses dt as a reaction/source interval. It is not the
transport timestep and must never be labeled as one. Reaction positions are
fixed during each source half, as in the original symmetric composition.
The atomic integrator still evolves actual local populations/momenta/energies
and predicted midpoint rates; this does not freeze reaction populations for the
whole run. Birth supports and Gauss/field changes retain the precise projection
contract. No physical-width coarsening parameter is available.

This trades source temporal resolution for fewer source evaluations; it does
not certify that a particular source interval is adequate. A symmetric split
has the usual formal second-order source/transport error in regular resolved
limits, but variable-space projections can add errors. Complete source-interval
refinement and local field/species comparisons are required before adoption.
The finest acoustic/plasma transport resolution is retained independently.

`last_transport` records every accepted internal dt, its current admissible
bound, and actual energy/thermal diagnostics. These are split-stage states,
with the first source half already applied and the second not yet applied;
they must not be presented as completed coupled states. Completed calls return
the full source interval and both physical source reservoirs. An invalid stage
or exceeded declared step budget rejects the full interval without shortening
the requested physical duration. This is an opt-in research composition.

`progress_every` optionally prints actual accepted transport time and step size
at that integer interval. It defaults to zero and changes no state, stopping
condition, timestep or physical source interval. The long argon driver uses100
to expose ongoing fine-resolution computation before a complete source return.
