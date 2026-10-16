# Conservative source-to-material projection

`project_reactions(overlap, state, ...)` returns a MaterialSpeciesSet and state
with actual changed particle amounts and newly occupied support. A block is
created for each connected nonempty interval sequence of each atomic species;
empty species and exterior intervals remain absent. No original boundary is
reused to discard product particles. No density threshold defines occupancy.

On the exact common overlap intervals, N_s=integral(n_s dx) and
E_s=integral(total translational energy density dx) are Gaussian-integrated.
All species share these partition edges, so locally conserved stoichiometric
charge has the same interval average before/after reactions. Original initial
charge is constant on each overlap interval; its Gauss field is unchanged up
to quadrature/source roundoff. No separate charge deposition or neutrality reset.

For each support, the consistent P1 mass matrix M uses m_s N_s/6[[2,1],[1,2]].
Nodal velocity solves Mv=integral(momentum_density phi dx). Partition of unity
preserves total vector momentum. Resolved cell kinetic energy K_s is then
integrated exactly and actual thermal energy is recovered as U_s=E_s-K_s.
This is a conservative moment representation change, not an energy-residual
correction: every cell retains its source-evolved total energy. No global budget
is inspected or used to add heat. Negative U rejects without clipping, floors,
redistribution or retries. Positivity of this projection is NOT guaranteed.

Momentum/energy projection stays in the atomic source's explicitly fixed
Galilean frame. The returned material set carries the same frame and restores
lab outputs. Exactly cold common translation therefore does not require a
subtraction of two large lab kinetic energies during thermal recovery.

This projection changes local velocity/current and can convert resolved bulk
structure into thermal energy; that is numerical mixing, not physical collision
heating. It must approach identity for a zero-duration source on representable
initial fields and converge when the actual source profiles are resolved. Source
generated discontinuities can be inadmissible in continuous nodal velocity;
such failures remain evidence that a richer representation is needed.

Photon, prescribed-heat and atomic-binding amounts are returned explicitly.
The per-block transport total excludes these reservoirs; callers must retain
them in the physical energy budget. No photon transport is implied.

Default `coalesce_ulps=0` retains every overlap edge. An explicit nonzero value
permits removing redundant internal cuts separated by at most that many ULPs
at the remote-domain coordinate scale. Endpoints of EVERY occupied support
are protected, including extremely thin species. All original quadrature
weights/moments are integrated into the merged intervals; no particle or energy
is discarded and no coordinate endpoint is moved. Returned records disclose
the actual tolerance and removed edge count. Local Gauss profiles can change
at roundoff scale; field-energy effects must be measured, never repaired.

This addresses floating duplicate cuts only. Distinct physically resolved cuts
still proliferate under repeated reacting/moving composition, and practical
coarsening/timestep control remains unresolved. The option is not a general
adaptive remap or a formation replacement.

Measured limitation: the source-only reacted fixture changes J by0.2237% during
the uncoalesced projection and0.3172% with8-ULP coalescing, despite roundoff global
momentum/energy differences. Removing a near-coincident cut changes the velocity
space and is not necessarily a roundoff change in current. The actual moving
argon attempt still collapses to~1e-25..1e-23s after its first step because
protected support boundaries and distinct internal cuts remain. Do not increase
this tolerance or remove support ends to conceal that failure. A separate
density/quadrature versus velocity space is assessed in
examples/material_velocity_resolution.md; it has not been implemented in motion.
