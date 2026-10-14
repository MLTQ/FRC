# Moving carrier-current scientific comparison

Reads completed reference and free trajectories; produces a scientific PNG/PDF
without re-running a simulation or adding an acceptance framework. The inertial
reference is compared with its independent amplitude law and the distinct ideal
flux-freezing prediction. Its prescribed motion has external mechanical work;
its total-fluid energy/virial residual is not an unforced conservation check.

Free-run magnetic/current/thermal energies and boundary evolution use saved
actual states. Toroidal magnetic field and poloidal current are reconstructed
through the same moving physical compatible spaces. Current is en*w from the
carrier-weighted Ampere solve, not division of a cellwise curl by a floor.
Cell colors are volume averages and current arrows are display samples.
The figure labels the isentropic toroidal sector and missing coil/formation scope.

The existing conditional Coulomb `regime` helper is evaluated at maximum-density,
maximum-drift and minimum-density quadrature samples of the displayed initial
and final states. The reference length is the body radius, not a resolved local
gradient length. Collision counts flag the eta=0 mathematical limitation rather
than silently installing a resistivity. No validity certificate is inferred.

Initial free-run B and J are compared against independently differentiated
compact cylindrical formulas, including particle-weighted drift error. A separate
80x80 physical integral gives reference magnetic and current kinetic energies.
These errors include the initial magnetic projection and current representation;
conservation does not normalize them away. Two free grids give a comparison,
not a demonstrated asymptotic convergence order for the nonlinear trajectory.

The independent initial acceleration also includes pressure, JxB and the full
mu*n*w*w stress divergence, with cylindrical connection terms. Complex steps
differentiate only the explicit analytic compact reference. The actual force is
recovered through the production consistent mass solve and compared in a
particle-weighted norm. Historical magnetic-only initial representations remain
in the output comparison, including their tradeoffs against the new physical
magnetic-plus-current projection.

The two spatial panels show the final Btheta and quantitative signed Jz, with
current-direction arrows and the initial boundary dashed. They do not rescale
an animation or imply a 3-D stability result. Initial states remain in the
independent comparisons and energy histories.
