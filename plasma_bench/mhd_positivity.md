# Admissible reconstruction helpers

`positive_scale` contracts a scalar reconstruction toward its positive cell
value. `positive_blend` finds a convex weight between an admissible thermal
pressure trace and a proposed total-pressure-derived trace. Both enforce a
small relative trace margin (1e-10 times the local positive cell value), with
roundoff slack. This is not a density/pressure floor in the conserved state.
No value in a cell's mass, momentum, energy or magnetic flux is changed.

`limiter_statistics` counts activations and records the minimum blend weight;
the dynamic solver accumulates accepted RK-stage diagnostics. An activation
is numerical limiting, not a physical heat or transport source.

These helpers bound reconstruction evaluation states only. They do not prove
positivity of the full multidimensional MHD time update, whose magnetic-energy
and discrete-divergence constraints require additional analysis. Invalid cell
states still raise. The accompanying method is independently derived, not a
claim to implement a published globally positivity-preserving CT theorem.
