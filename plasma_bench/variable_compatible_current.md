# Variable-density compatible current energy and induction

`VariableCompatibleCurrent` supplies the existing original-equation induction
solver with the same physical M, H and Faraday C, but a variable-density inertia
I=P.T W [mu/(n e²)] P. Here n=R n_cell uses actual conservative cell averages;
P samples the compatible electric/current fields. Positive quadrature weights
and positive n make I symmetric positive definite. Its tensor matrix-free action
avoids storing a large three-dimensional point map. `weak_vector` is P.T W.
Physical electron and ion masses are unchanged. No frozen density is substituted.

Ampere remains M J=−C.T H deltaB. The immutable affine field retains its full
uniform bias in magnetic energy; current acts on its independent increment.
The induction solver can use I directly without changing its original equations
or tolerance. An optional existing geometry avoids rebuilding its fixed factors.

For actual supplied n_t=R(rho_t/(mi+me)), the explicit Ohm forcing is
F_n=−P.T W [mu J n_t/(e² n²)]. The local current-energy derivative is

    K_J,t = integral_cell mu/(n e²) [J dot J_t − |J|² n_t/(2n)] dV / V.

`current_energy_rate` exposes both contributions separately. The global
`density_energy_gradient` is R.T W[−mu |J|²/(2 e² n²)] and must be paired with
the actual continuity rate. It is not simply a nodal or cell-center potential.
`magnetic_energy_rate` retains the local B0 dot deltaB_t cross term even though
its global integral vanishes at these fixed-flux boundaries.

An induction solve with force F_other+F_n obeys

    d(W_B+K_J)/dt = −J.T F_other + integral mu |J|² n_t/(2 e² n²) dV.

The last term is an exchange that a full moving energy/pressure update must
account for, not external heating. This component does not evolve n, velocity,
species entropy or total energy and supplies no missing pressure-work closure.
All boundaries remain the bounded compatible-field fixture, with no vacuum,
diffusion, density floor, adjustable inertia or energy reset.
Unrepresentable positive inertia (for example density so tiny that n e²
underflows) rejects explicitly; it is never converted to zero or infinite mass.
