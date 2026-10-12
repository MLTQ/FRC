# Complete workflow with componentwise drift stress

`run` temporarily dispatches the scheduled signed workflow to
ComponentStressTransport. All original inputs, equations, solver checks,
independent q16/20 volume/face audits and early raw records remain. The complete
explicit source closure includes this wrapper and its selected implementation.
Serial bindings restore on success or failure; do not run concurrently in one
process. Fresh paths only. Complete old/new comparisons and fresh full replays
must establish numerical equivalence; source similarity alone is insufficient.

No unit tests, new diffusion, floors or physical-model changes. Measured memory
and numerical fidelity are separate outcomes. This is an instantaneous
verification fixture, not a time trajectory or physical FRC experiment.
