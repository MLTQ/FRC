# Magnetic/vacuum work subsystem

`MagneticStepper` advances induced nodal flux u and fluid momentum using implicit
midpoint. The exterior equation is K u=0, not a low-density or low-conductivity
fluid. Let g=(2pi/mu0)K u, let C(psi) interpolate v·grad(psi) to plasma nodes,
and let Mv contain actual fluid cell masses. The paired equations are
u_dot=-C v-psi_external_dot-R K u and Mv v_dot=C^T g. Here
R=eta*r/(mu0*plasma_node_area) in plasma and zero outside. Vacuum rows are
replaced by their elliptic constraint. No plasma mass exists there.

The nonlinear midpoint is solved by iterating its sparse, velocity-eliminated
linear system. Nonconvergence raises; no partially updated state is returned.
The gas energy increment is the actual kinetic work plus positive Ohmic heat,
not a repair to a failed total-MHD-energy update. Heat uses nonnegative cell
shares; the complete field/gas work identity is independently diagnosed.

`MagneticUpdate` returns fields, gas, integrated source/Ohmic/Lorentz work,
exchange residual, nonlinear residual and iteration count. Source work is
-g_mid·delta(psi_external). The energy is gas plus induced self-field energy;
external-field self/mutual energy is not mislabeled as part of that budget.
See docs/COUPLED_NUMERICS.md for derivation, assumptions and limitations.

Adjoint pairing proves discrete work cancellation, not spatial force accuracy,
shock fidelity or free-boundary validity. These require independent checks.

With multipole boundaries, the same K includes unbounded exterior magnetic
energy. Algebraic vacuum rows also include the modal amplitudes; C, R and source
forcing have zero modal rows. Prescribed coil/bias fields are sampled on physical
nodes, never projected into the decaying induced-field basis. The paired work
identity therefore remains unchanged, including exterior energy beyond the grid.
