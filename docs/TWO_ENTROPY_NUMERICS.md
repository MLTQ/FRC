# Joint physical energy, momentum and species entropy transport

The common radial-face prototype corrected near-axis forces but retained a
spatial convective energy defect. A force-only correction was already shown to
violate either original force accuracy or uniform-translation momentum. This
separate component instead derives a joint flux for all transported variables.
It does not change the frozen radial-face study or replace the plasma solver.

## New state and physical equations

The state is the point value U=(n,mx,my,mz,C_e,C_i), where m=n v and C_a=n s_a.
These are not cell means or the previous independently nodal specific entropies.
For smooth ideal species, n, all Cartesian momenta and each C_a obey conservative
advection, with total pressure in the momentum equation. In SI units,

`h_a = e exp(s_a) (n/n0)^(gamma_a−1)/(gamma_a−1)`,
`p_a = (gamma_a−1)n h_a`,
`H = M |m|²/(2n) + n(h_e+h_i)`, with M=mi+me.

The physical log entropy per particle differs from s_a only by a fixed positive
species factor and additive reference constant. Conserving each C_a therefore
conserves each physical species entropy in this smooth adiabatic model.
No energy account is used to reset or infer either temperature.

## Two-state derivation

The physical H gradient is

`w_n = −M|v|²/2 + sum_a (gamma_a−s_a)h_a`,
`w_m = M v`, and `w_Ca = h_a`.

The corresponding potential for direction d is psi=p v_d. Let bars be
arithmetic means, a superscript log be the positive logarithmic mean, and
brackets be right-minus-left differences. Choose

`Gamma = n^log bar(v_d)`,
`F_m = Gamma bar(v) + bar(p)/M e_d`,
`F_Ca = bar(v_d)[n^log bar(s_a) + (gamma_a−1)(bar(n)−n^log)
                  + n^log(bar(h_a)/h_a^log−1)]`.

The kinetic terms cancel by the arithmetic velocity average. The exact
thermodynamic jump relation

`[s_a] = [h_a]/h_a^log − (gamma_a−1)[n]/n^log`

then gives `[w] dot F = [psi]`, separately accounting for both species before
summing. The symmetric local physical-energy flux is
`G = bar(w) dot F − bar(psi)`.

This derivation is for our state with two entropy densities. It is not a claim
that an existing paper already supplies a verified two-temperature FRC solver.
The kinetic-preserving momentum form follows the established construction in
[Chandrashekar (2013)](https://arxiv.org/abs/1209.4994). Centered high-order
combinations of two-state fluxes are an established framework; see
[Klein and Öffner](https://arxiv.org/abs/2211.01171). Our new formula and its
physical two-species identities are checked independently in code.

This flux does not exactly preserve constant specific entropy pointwise when
density varies. Its correction to F_C differs from Gamma times that constant;
local entropy variations are truncation error even though the global entropy
sums are conserved. The isentropic wave therefore requires separate local
entropy convergence measurements, retained in TWO_ENTROPY_VERIFICATION.

## Spatial and temporal conservation are distinct

On periodic uniform point nodes, the skew derivative D has zero row sum.
The update `U_t,i = −2 sum_j D_ij F(U_i,U_j)` conserves all six state sums.
The jump identity proves the local physical-energy equation
`w_i dot U_t,i = −2 sum_j D_ij G(U_i,U_j)`. Thus physical H, not an independent
diagnostic account, is conservative in the declared nodal dx integration norm.
Orders2/4/6 use their ordinary centered derivative coefficients.
Axial integrals are per unit cross-sectional area: H has units J/m³ and its
reported integral has units J/m². No cylindrical volume is inferred here.

The theorem concerns that quadrature. Independent continuous physical energy
is obtained from a declared real trigonometric interpolant of U, recomputing
n, velocities, both temperatures/pressures and H on denser nodes. Its difference
from nodal H and its own drift are retained. A quadrature identity is not an
excuse to omit physical integration accuracy.

RK4 is not exactly energy preserving. Fixed-grid step-halving tests measure its
remaining energy error; no post-step correction is applied. Every stage and
final state must have positive representable physical density and temperatures.
Failure leaves the accepted state unchanged. There is no artificial diffusion,
limiter, mass/temperature floor or shock entropy-production theorem. The state
is positive-density bulk gas, with no exterior atmosphere introduced.

## Required cylindrical and electromagnetic work

This is an axial transport verification, not the completed moving plasma model.
The next derivation must use the actual radial field, cylindrical volume,
angular momentum n r vtheta, regular-axis space and impermeable-wall response.
Pressure, centrifugal and geometric torque work must preserve both physical
energy and conservative angular momentum while matching the original Euler
forces. A point-state axial result cannot be copied into the old cell/face
state or assigned its energy theorem without proving those maps.

Verify the same independent regular polynomial force, full two-dimensional
RadialWorkReference rates, first-annulus/max errors and radial Bessel modes.
Then repeat unforced multidimensional dynamics and the complete current-energy,
Maxwell/drift, species-enthalpy/Poynting coupling with actual continuity n_t.
The old variable-density current operator expects its declared density space;
any new bridge must preserve its physical energy/continuity contract explicitly.

Heated free-edge/exact-vacuum treatment, physical collisions/transport/ionization,
shocks, real boundaries and prescribed coils, converged formation and experiment
controls remain requirements. This component does not demonstrate confinement
or three-dimensional stability. Track the joint integration as frc-7mm.1.
