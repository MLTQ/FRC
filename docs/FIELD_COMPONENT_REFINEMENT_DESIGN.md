# Independent field-direction refinement for current-component accuracy

Continue from verified checkpoint09b9363, classified as progress. No unit tests
are added or run. Keep the full original FRC sandbox scope and12 acceptance
beads open. Earlier source archives and the user's poster remain untouched.

The fixed independent radial16,N128 contact-mode reference and F64/256 fluid
representation remain unchanged. Baseline fields are B256/256, epsilon1e-4,
induction tolerance1e-10. Independently double the axial field resolution to
B256/512, then radial resolution to B512/256 where resources permit. Every
case retains physical masses, actual density, all current/species terms,
original residual definitions/bounds/budgets and independent q16/20 audits.

Source-derived storage estimates precede execution. Run only one large
geometry/solver process at a time, retain exact inputs and progress, and guard
disk headroom. Guard interruption is a resource event, not a solver rejection.
Never turn a failed derivative into accepted data or weaken bounds to force a
passing refinement. Original/replay builds use separate processes.

The original axial refinement was resource-interrupted during initialization,
before any signed input or derivative was saved. Preserve this partial run;
do not replay resource exhaustion merely to obtain a profile. Its own peak RSS
reached 11.23 GiB, while free disk crossed the 256 MiB guard. Neither a numerical
rejection nor a field-accuracy result exists for this attempt. Native work did
not yield promptly to SIGINT; approved termination was collected as exit 143.

The opt-in tensor geometry evaluates each electric mass inverse using separate
radial and axial LU factors. For M = Gr kron Gz and row-major vectorization,
M vec(X) = vec(Gr X Gz^T). This preserves the original basis, assembled mass,
Hessian, curl, quadrature and boundaries; all original physical residuals use
the original assembled matrices. Finite-precision evaluation differs, so full
original/new physical comparisons and strict replays are required. The default
geometry and all historical Python sources remain unchanged.

Complete coarse signed runs observe both actual base and strong-wall operators;
all operator fingerprints must agree and each full run must reproduce its own
archived pilot. Fine-grid signed runs compare all physical components with the
original baseline, retaining individual current errors and actual differences.
These are numerical-equivalence controls, not evidence that spatial error fell.

A complete fine replay observes entry/exit cumulative own-process peak RSS in
current transport, species pressure and metric construction. These are read-only
observations with a separate source archive. Require every scientific array and
metric to reproduce the unobserved baseline. Cumulative high-water intervals
are neither live/exclusive allocations nor an additive memory budget. The new
guard enforces 8 GiB own peak RSS and 256 MiB free disk; if SIGINT does not yield
within five seconds, only its own worker exits 130. Larger-grid feasibility is
still unproved; target measured costs before another full refinement.

Completed, independently computed replay NPZ files may share storage only after
their verification record and both complete byte hashes are checked. Preserve
all paths/bytes and record each atomic hard-link replacement. Shared artifacts
are immutable: neither pathname may subsequently be modified in place.

The existing complete F32/256,B256/256 signed case is independently reconstructed
and fully replayed for component/strong-curl norms. It provides a fixed-field
fluid-resolution comparison to the accepted F64 baseline before allocating a
new high-cost radial field solve. Original saved rates and hashes remain intact.

Accepted cases require complete independent component norms, actual Bt strong
curl/current comparison, quadrature controls and complete fresh signed replay.
Same physical reference and component normalizations permit directional
comparison despite different meshes. Strong curl excludes distributional
face jumps; local conservation and complete time trajectories remain separate
requirements. Rejected cases require exact preserved inputs and full failure
replay before diagnosing arithmetic. Refinement cannot be inferred from a
small total-current norm alone.

The parent equations, periodic filled-argon domain, impermeable wall and
constant prescribed bias are unchanged. No numerical diffusion, filter,
atmosphere, density floor, artificial mass or energy repair is introduced.
These fixtures do not select lab conditions, prove formation, establish3-D
stability or show confinement. Exact moving vacuum, transport/ionization,
physical coil/end boundaries and experimental UI remain open.
