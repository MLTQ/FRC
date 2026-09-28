# Dynamic MHD verification
Exact annular geometry and uniform stationary balance cover the coordinate axis.
A one-dimensional reduction checks that UCT retains the full Rusanov magnetic
dissipation. Independent cylindrical fast-wave, translated 2-D weak-flux and
exact Sod fixtures measure refinement of motion, induction and shocks.
A weak-field Bessel decay tests physical resistivity, positive thermal heating
and energy partition without adding energy. Closed-boundary mass/energy and
actual face divergence are checked during dynamic runs. Invalid density,
pressure, unsupported radial boundaries and excessive time steps must fail
explicitly without modifying the accepted state.
Unsupported-wall fixtures reconstruct total energy at positive pressure so the
boundary contract, rather than a separately invalid energy state, is tested.
At a fixed spatial grid, step refinement against a much finer time reference
independently verifies second-order RK2 time accuracy without confusing it with
the first-order spatial error.
The isolated magnetic acceleration is checked against analytic J×B independently
of gas pressure. A static total-pressure equilibrium starts force balanced;
subsequent artificial flow caused by numerical diffusion must decrease with
refinement. This does not claim exact equilibrium preservation or FRC balance.
Magnetic-force convergence uses the physical annular volume norm and includes
the axis cell. A grid-index-trimmed unweighted norm would change its physical
region with refinement. The separate report retains first-axis-cell absolute
acceleration error; it decreases only linearly for this first-order method.
