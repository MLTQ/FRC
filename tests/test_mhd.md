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
spatial error. Smooth early-time refinement uses both spatial orders, while
the report separately retains the less regular quarter-period limiter-switch
measurement; its temporal convergence need not have clean fourfold ratios.
The isolated magnetic acceleration is checked against analytic J×B independently
of gas pressure. A static total-pressure equilibrium starts force balanced.
First-order artificial flow must decrease with refinement, while the reconstructed
polynomial fixture must stay near roundoff. This does not claim general FRC balance.
Magnetic-force convergence uses the physical annular volume norm and includes
the axis cell. A grid-index-trimmed unweighted norm would change its physical
region with refinement. The separate report retains first-axis-cell absolute
acceleration error: linear for the reference method and cubic for the new
axis-regular reconstruction on this polynomial fixture. Polynomial annular
averages independently test the axis parity and moment reconstruction.
Both first- and second-order UCT reductions must retain full 1-D dissipation.
The second-order constant-total-pressure straight-field fixture is preserved
to roundoff; the first-order baseline still has measurable decreasing drift.
