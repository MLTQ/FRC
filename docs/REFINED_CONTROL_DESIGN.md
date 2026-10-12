# Refined-grid amplitude and residual controls

Continue from cd66ec5, verified progress: all preceding workers terminal,
all source/evidence preserved. No unit tests. The full FRC goal and all twelve
acceptance beads remain open.

Use the accepted component-stress physical workflow unchanged, with F64/256,
fixed independent radial16/N128 reference and separate B256/320, B320/256 grids.
Explicitly vary centered perturbation epsilon from the existing1e-4 baseline
to1e-3 and1e-5, retaining both signed full nonlinear rates and q16/20 physical
and face audits. Compare actual component differences as well as their errors;
the full-current norm hides weaker components. Reconstruct the actual magnetic
rate's strong curl independently, without dropping face-distribution caveats.

At fixed epsilon1e-4, request a tighter1e-11 original-equation tolerance on
each refined grid. Accepted rates must meet it. Strict numerical failures keep
their exact inputs and full rejection metadata and require fresh complete
replay; no partial iterate is accepted and no bound is weakened. Existing
raw-curl/addition-loss evidence is not presumed repaired by spatial refinement.
Any successful tighter-bound control also requires full physical component
and actual magnetic-curl comparison against its same-grid original-bound
baseline, with complete independent audit replay. Rejections have no physical
error curve and remain explicitly recorded.

Every new calculation runs in a separate process, under unchanged8GiB own-RSS
and256MiB free-disk guards with final-extrema checks. Do not run another large
solver or independent diagnostic concurrently. All completed physics and
rejections require fresh whole-protocol replay; resource interruptions require
terminal process collection and preserved partial evidence before fresh retries.

Freeze new control sources and companions before calculations; do not edit
already used Python files. Preserve every previous artifact and the user's
poster. Byte-identical completed verified records may share immutable storage,
with full hash/verification manifests; never overwrite a shared path.

These filled-periodic wall fixtures measure consistency and sensitivity, not
time accuracy, formation,3-D stability or confinement. No added diffusion,
filters, massive/conducting atmosphere, artificial inertia, density floor or
heat/energy repair. Continue combined/further spatial controls, full signed
trajectories and local budgets, then original heated exact-vacuum, physical
transport/ionization/validity, prescribed ends/coils, formation and UI gates.
