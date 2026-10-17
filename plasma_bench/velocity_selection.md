# Initial velocity knots relative to requested targets

The nearest-density-face map can assign neighboring nominal cosine targets to
different faces within the same tiny density cluster, creating a mechanical
interval many orders smaller than requested. This explicit research alternative
first computes that usual map, then resolves target collisions before any new
velocity state is created. A candidate represents its closest associated target.
If two candidates are closer than a declared fraction of their associated target
spacing, retain the better-positioned one; physical endpoints always remain.
Repeat until the selected mesh meets that target-relative separation.

The fraction is dimensionless, defaults to zero (the historical map), and must
lie in[0,1]. The experiment uses0.5. Actual velocity-cell counts can be lower
than the requested count and must be reported; this is a change of the initial
approximation space, not a guarantee of local accuracy. No existing independent
velocity knot passes through this selection: retained meshes bypass it. Every
fine density cell and physical endpoint is unchanged, with the same N/U/mass.
There is no minimum physical density width, force modification or timestep fix.

Full charged runs and analytic expansion must evaluate this choice. Convergence
and local thermal/front errors remain requirements; avoiding tiny artificial
mechanical intervals alone does not validate a reacting plasma model.
