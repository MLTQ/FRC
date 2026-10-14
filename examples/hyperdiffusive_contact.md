# Moving contact with explicit numerical damping

Runs the complete original negative 1e-4 continuum perturbation for 100 us,
using the existing trajectory/8 GiB own-RSS guard and saved-state protocol.
Default F8/64, B32/64, 60 steps directly matches the earlier failing case.
No automatic replay, retries, relaxed bounds or unit tests are added.

The only equation change is HyperdiffusiveWallPlasma with explicit numerical
electron heat fraction 1. The fixed coefficient is
eta6=mu0 strength (c_fast+D_H/h) h^5/pi^6, strength=0.25 by default;
D_H=|(mi-me)/(mi+me) B0|/(mu0 e n0), c_fast=sqrt(cs^2+vA^2), and h is the
smaller mean radial/axial field spacing. This is a prospective mesh-based
choice, not a fit to an observed growth rate. The radial field grid is uniform
in r and nonuniform in x=r^2/2. The running source metadata incorrectly calls
it nonuniform in r; the stored coefficient uses the correct R/nr spacing. Its coefficient and all inputs are
saved in numerical_diffusion.json; eta6 vanishes with refinement. It is not
physical resistivity, a collision closure or a guarantee of adequate damping.

Use the existing trajectory_density_prefix observer for all accepted endpoints
against the unchanged independent continuum mode. Error includes numerical
damping as well as nonlinear, temporal and other spatial effects. A successful
run does not establish convergence or formation. Gauss high-frequency damping
limits remain explicit in the solver companion.

Temporary factory bindings reuse the original driver in one process and are
restored in finally blocks. Fresh physical/control source archives include this
module and its new model/solver. Archived original modules remain unchanged.
