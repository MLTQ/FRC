# Mechanical refinement without deleting material cuts

`refine_velocity_faces` enriches an existing nested P1 mechanical mesh using
the same precise fine density faces. For each old velocity element, select the
existing density face closest to its midpoint (lower face breaks an exact tie).
Add it only within the middle half of that element. Every old knot remains;
the new element is at least one quarter of its parent length. No density face,
physical endpoint, population or heat is deleted or coarsened.

Repeat the declared number of rounds. If a complete proposed round exceeds
the declared mechanical cell capacity, skip that whole round and stop. This
keeps selection independent of traversal order and records the capacity limit.
Existing meshes above the capacity are retained unchanged. Missing central
density faces can prevent subdivision; actual counts, additions and completed
rounds are reported, not inferred from the requested level.

This is an explicit approximation/resolution control, not a physical transport
coefficient or an error estimator. It does not promise a requested accuracy.
Keeping density discontinuities while allowing a richer velocity basis targets
the measured trace-ion heating: in the actual 10–15 ns transport stage, refining
only Ar8 from 32 to 44 velocity elements reduced maximum temperature from
7.0414 to 0.3909 eV. Both full20ns all-species and heavy-only enriched10ns-source
runs now complete. The original5ns-source attempt fails on cold-cell arithmetic
and is preserved separately. Both corrected full5ns-source cases now complete,
but final Ar9 reaches370.17eV despite small total-energy residuals. Earlier Ar8
improvement does not establish source-time thermal convergence. The actual final
transport-stage refinement in outputs/reaction_final_trace_transport measures
this new local heating against independent adiabatic compression.

Source return constructs its ordinary inherited/newborn mesh first and then
enriches it, retaining the existing support-extension setting. It solves the
selected physical moment/energy transfer in that enriched space. Subsequent
motion uses its actual acoustic bound; no enlarged step or hidden retry.
Defaults request zero rounds and preserve previous behavior. Refinement can
increase cost, change numerical diffusion and still leave unresolved regions.

Inputs: ordered precise density and velocity faces, nonnegative integer rounds
and positive integer cell capacity. Output: old plus new velocity faces and
the explicit refinement record. It does not modify its input arrays.
