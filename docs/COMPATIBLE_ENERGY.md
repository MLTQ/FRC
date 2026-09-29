# Compatible staggered energy and magnetic-vacuum work

This is a selectable extension of the [conservative fixed-chamber reference](CONSERVATIVE_COUPLING.md).
It addresses a mismatch between magnetic energy inside the fluid and the
variational energy used by the exterior. It does not introduce an atmosphere,
change the independent shock fixture, repair gas temperature, or establish
physical formation. The historical method remains the default comparison while
the new controls undergo spatial, temporal and driven-source verification.

## Positive fluid magnetic energy

For a cell with radial midpoint r_c and face radii r_-, r_+, define

```
M_cell = [ (r_- Br_-² + r_+ Br_+²)/(2 r_c)
           + (Bz_bottom² + Bz_top²)/2 ] / (2 mu0).
M_fluid = sum(V_cell M_cell) = (pi/mu0) psi_total^T K_fluid psi_total.
```

All weights are nonnegative. The axis radial-face contribution has zero weight.
K_fluid contains precisely the fluid portions of the variational magnetic edges;
interface dual edges have half measure. Squaring interpolated cell B instead
loses short-wavelength magnetic energy. Both conventions converge for smooth
fields, but they cannot be interchanged after advancing a total-energy state.
Construction, pressure recovery, admissibility and diagnostics select the same
`energy_quadrature`. Face CT fields and cell vector fields for momentum are still
distinct numerical representations of the same continuum field.

The pressure condition E - momentum²/(2 rho) - M_cell > 0 remains a concave
constraint for rho>0 because face fields depend linearly on psi. Checked Euler
proposals therefore admit SSP convex combinations. This is a conditional
admissibility argument, not a general baseline positivity or entropy theorem.

## Natural vacuum traces and actual electric work

Let u be full-space induced flux, including harmonic vacuum completion and
decaying multipoles, K its symmetric stiffness, and psi the total fluid flux.
The response-energy derivative after subtracting prescribed source power is

```
dE_response/dt - P_source
  = -P_out + (2 pi/mu0) [restrict(K u) - K_fluid psi] dot dpsi/dt.
G = restrict(K u) - K_fluid psi; dpsi/dt = -r E_theta.
```

At noncorner boundary nodes use Bz_wall=-G/dz, Br_bottom=-G/dr and Br_top=G/dr.
These are the natural discrete boundary loads, not a redistributed global
residual. At an outer/end corner, half of each boundary length contributes to
one weak equation. Independent second-order exterior derivative estimates are
projected onto it by minimizing the boundary-length-weighted squared change of
the two traces. Harmonic polynomial tests check their continuum convergence,
including corners and the regular axis.

Impermeable-wall Poynting flux integrates the nodal products of the **selected**
CT electric field and these traces. It does not multiply separately averaged
E and B. The two-node products carry their radial/axial dual measures, so their
boundary sum equals `(2 pi/mu0) G_boundary dot dpsi/dt`. The resulting face flux
enters the ordinary conservative divergence. Maxwell stress uses corresponding
face-averaged magnetic traces. No state energy is subsequently modified.

First- and second-order fluxes at a stage share the same initial G. Local
fallback reassembles Poynting flux after selecting shared electric corners;
every flagged boundary cell then has the checked first-order face work exactly.
The conditional neighbor-recheck argument remains applicable.

Interior G is **not** moved to the boundary. With constant uniform bias it
vanishes to roundoff, giving a semidiscrete response identity. Analytic coil
fields generally have finite discrete K_fluid residual despite being continuum
vacuum fields. The remaining interior source-discretization error is measured
under mesh refinement. No general exact source-energy theorem is claimed.

## Time integration and a separate finite-step error

The optional `time_order=3` uses the three-stage SSP method, with RHS evaluation
times t, t+dt, t+dt/2 and budget weights 1/6, 1/6, 2/3. Each stage advances total
fluid flux, subtracts the analytic external field at its actual stage time and
completes the induced vacuum field. Every Euler proposal and final state must
validate before state, counters or work commit. The standard SSP construction
is described by [Gottlieb, Shu and Tadmor](https://www.math.umd.edu/~tadmor/pub/linear-stability/Gottlieb-Shu-Tadmor.SIREV-01.pdf).
Independent nonautonomous manufactured equations verify second/third temporal
order and exact budget quadrature with actual changing source fields.

SSP RK2 and RK3 do not generally preserve quadratic invariants exactly. For
constant bias, write response energy as a linear function plus psi^T H psi,
where H=(pi/mu0)(extension^T K extension-K_fluid). Direct expansion gives

```
Delta E_response
 = dt sum_i b_i spatial_defect_i
 + dt² sum_ij (b_i b_j - b_i a_ij - b_j a_ji) (k_i^T H k_j).
```

Here k_i are total-flux stage rates. This is the usual quadratic-invariant RK
coefficient condition; see [Sanz-Serna's analysis](https://arxiv.org/abs/1503.04021).
The audit independently completes each rate through vacuum and contracts this
form; it never uses the observed energy change to construct its prediction.
Agreement separates finite-step error from a spatial work mismatch. It is not
an energy correction or a claim of time-converged boundary dynamics. Local
limiter choices themselves depend on timestep, so isolated CFL halving need
not produce monotonic error at an unresolved sharp boundary.

## Verification contract and remaining scope

Independent tests cover uniform energy, analytic polynomial energy and traces,
the instantaneous no-spurious-heating identity at J=0, evolved-shock work,
nonautonomous stage order and atomic rejection. Recorded studies retain the
original shock, source pulse, smooth analytic solutions and method ablations.
They report source-normalized response error, fluid-flux conservation, physical
and numerical heat, positivity, CT, vacuum residual and actual limiter stages.
The [measured report](COMPATIBLE_ENERGY_VERIFICATION.md) distinguishes improvement
from unresolved time/mesh sensitivity. Free surfaces, actual gas ends, transport,
Hall/kinetic validity and converged formation remain separate requirements.
