# Guard complete detached-cache refinement

Uses the unchanged fine-workflow settings and resource guard: F64/256 fluid,
specified field nr/nz, radial16,N128 reference, epsilon1e-4, tolerance1e-10,
own-worker8GiB peak RSS,256MiB free disk and5s bounded self-termination after
SIGINT. Only the detached-cache driver is selected. Control archives explicitly
record that dispatch and its full source closure.

Fresh paths are mandatory. --original verifies complete same-source scientific
replay; algorithm comparisons use the whole-workflow comparison helper.
Temporary bindings restore on all exits. Resource events do not change numerical
acceptance, authorize a failed rate or establish refinement feasibility.
No unit tests or old source modifications.
