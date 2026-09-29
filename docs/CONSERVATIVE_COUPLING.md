# Conservative MHD with a magnetic-only exterior

`ConservativeVacuumMHD` advances mass, momentum, total fluid-cell MHD energy and
magnetic flux together. It addresses the split reference's inconsistent fluid
and magnetic shock diffusion. The original split solver remains available for
its exact work pairing and comparison; neither solver yet provides a moving
free plasma surface or validated FRC formation.

The equations below describe the historical/default cell-energy, Riemann
boundary-work and RK2 choices. The selectable staggered-energy, natural-boundary
and RK3 extension has its own [derivation and limitations](COMPATIBLE_ENERGY.md)
and [measured study](COMPATIBLE_ENERGY_VERIFICATION.md). Existing study drivers
now select their historical controls explicitly.

## Equations, geometry and state

Inside the fixed gas-filled r-z cylinder, the equations are compressible
resistive MHD with no toroidal velocity or magnetic field:

```
d rho/dt + div(rho v) = 0
d(rho v)/dt + div[rho v v + (p+B²/(2mu0)) I - B B/mu0] = 0
d E/dt + div[(E+p+B²/(2mu0))v - (v·B)B/mu0 + eta J×B/mu0] = 0
d psi/dt = -r E_theta; E_theta = vr Bz - vz Br + eta J_theta
E = p/(gamma-1) + rho v²/2 + B²/(2mu0)
```

Here J is electric current density in A/m². Consequently the resistive energy
flux is **eta J×B/mu0**, since E_res=eta J. Cylindrical finite-volume radial
momentum includes the integrated total-pressure metric source. Flux convention
remains Br=-psi_z/r, Bz=psi_r/r, with regular zero psi on the axis. Shared face
flux differences give compatible divergence at floating-point accuracy.

`ConservativeState` stores cell total MHD energy and the full induced flux. This
differs from the split `CoupledState` gas-only energy. Conversion is explicit and
allowed only for an unadvanced initial fixture, not as an evolved-state repair.
No magnetic energy is double counted in either formulation's diagnostics.

External coils and bias are still analytic `VacuumField` sources. Induced flux
outside the gas obeys the existing symmetric elliptic zero-current constraint,
with the verified decaying multipoles for infinity (or explicit fixed comparison).
There are no exterior fluid mass, pressure or conductivity unknowns. Plasma-node
total flux evolves by CT. At each RK stage its exact external-flux increment is
subtracted and the remaining induced field is completed through the vacuum solve.

## Consistent interior numerical fluxes

Both fluid and magnetic fluxes use the same limited second-order primitive and
staggered-field reconstruction. Fluid fluxes are Rusanov; corner electric fields
use the four-state local Lax–Friedrichs upwind constrained-transport expression.
Its full one-dimensional limit matches the fluid Rusanov induction flux. A
separately selectable piecewise-constant reference uses the same construction.
No independently tuned volume shock-diffusion sensor is added in this solver.

Consistency between fluid Riemann fluxes and constrained induction is the central
UCT principle in [Londrillo and Del Zanna (2004)](https://arxiv.org/abs/astro-ph/0310183)
and [Mignone and Del Zanna (2020)](https://arxiv.org/abs/2004.10542).
This code's cylindrical metric, reflecting chamber, magnetic vacuum and
interface treatment are independent adaptations. Published Cartesian methods'
theorems are not automatically inherited.

In strictly perpendicular one-dimensional ideal flow, B/rho is advected with the
material. Shared mass/induction diffusion preserves a spatially constant ratio
in the interior to roundoff apart from actual two-dimensional boundary influence.
That particular invariant is not asserted for arbitrary curved two-dimensional B.

## Magnetic traction and energy at fixed gas boundaries

Gas uses reflecting slip traces: zero normal flow, freely reconstructed tangential
velocity. The boundary is magnetically transparent. Normal B is the shared CT
face field. Tangential B for Maxwell stress and Poynting transport is evaluated
from second-order **vacuum-side** normal derivatives of total psi. Independent
harmonic-polynomial tests check those traces, including the axis end faces.

Reflecting interior tangential B instead would omit the force/work of interface
currents. Even with exact fluid-domain energy conservation, the resulting whole-
space response error need not converge. Including exterior traction is necessary.
It also exposes strong, initially unresolved interface currents in the original
standing-shock fixture; the interior exact shock is not a globally equilibrated
plasma/vacuum configuration.

An explicit numerical interface diffusivity regularizes these unresolved surface
currents at **plasma boundary nodes only**:

```
a = |v| + sqrt(gamma p/rho + B²/(mu0 rho)); h = min(dr,dz)
eta_interface = mu0 a h / 4
J_node = (K u)_node / (mu0 plasma_node_area)
E_interface = eta_interface J_node
```

Neighboring positive cell averages supply a at nodes. A half-cell boundary dual
area makes this approximately the LLF electric dissipation a jump(B_t)/2. Its
sign comes from the positive stiffness current, not an indefinite difference of
one-sided derivative estimates. Vacuum, modal and interior plasma nodes have
exactly zero interface eta. The coefficient vanishes as h; this is numerical
diffusion and **not a material conductivity or reconnection model**.

The same electric field enters induction and Poynting flux on every adjacent
fluid face, including tangential transport along the boundary layer. Restricting
this term to outward energy flow would omit internal redistribution.
Independent positive numerical work is integrated as sum(eta_interface J² Vdual)
and reported separately from physical Ohmic heat. Neither is added as a separate
source to total energy: doing so would double count resistive energy transfer.
Finite-grid internal-energy change need not equal these heat diagnostics exactly.
The comparison switch `interface_diffusion=False` is retained, including failures.
Explicit wave plus physical/numerical diffusion bounds limit the time step.

## Complete conservative update admissibility

Time integration is SSP RK2. Positive reconstructed traces alone do not guarantee
positive updated pressure. Each full Euler proposal is checked before use.
If a second-order proposal is inadmissible and the first-order proposal is
admissible, the default local limiter replaces high-order face fluxes and corner
electric fields incident to every troubled cell, and that cell's metric source.
The selected shared fluxes are assembled again, including the actual boundary
energy rate. No cell value is clipped or independently reset.

The mixed update can trouble neighbors. Every new invalid cell is added and its
incident entities replaced, until the complete proposal is admissible. Already
flagged cells have exactly the checked low-order update. Every unsuccessful
pass adds at least one cell, giving finite termination conditional on a valid
baseline; there is no iteration cap or hidden global fallback. A shared corner
can affect four cells, so the reported affected region includes diagonal
neighbors. Cells beyond that region retain the original high-order derivative.
For translationally invariant radial flow, face/corner choices coincide and
preserve initially constant Bz/rho. General 2-D flow requires separate checks.

This is a local binary order reduction, not second-order evolution in troubled
regions. It follows the conservative a posteriori flux replacement principle
described, for different DG/subcell discretizations, by
[Vilar and Abgrall](https://doi.org/10.1137/22M1542696). Our cylindrical FV/CT,
metric and vacuum adaptation is separately implemented and checked; it does
not inherit their algorithm or its theorems. MHD positivity depends on the
particular discretization and divergence conditions; see the independent
analysis of [Wu](https://arxiv.org/abs/1802.02278).

For rho>0, internal energy
`E - momentum²/(2 rho) - B²/(2mu0)` is jointly concave in (rho,momentum,E,B).
Cell B is linear in corner psi. The retained comparison setting
`flux_limiter="global"` uses the interval of admissible convex weights from the
low-order proposal. Bisection finds its limit and uses 0.9 of that weight to
stay inside without an absolute pressure floor.
The two accepted Euler proposals then have an admissible SSP convex average.
If the low-order proposal is also nonphysical, the step raises; no general
positivity theorem for this plasma/vacuum interface is claimed.

Diagnostics record limited Euler-stage count, minimum high-order weight,
troubled/affected cell-stage counts, maximum affected fraction, replacement
iterations and per-cell cumulative masks. The global comparison can reduce
accuracy throughout the domain because of one boundary cell; local replacement
removes that direct loss in distant cells. It does not guarantee monotonic
whole-space response-energy error or eliminate indirect boundary influence.
Disabling the limiter retains the original cold-boundary failure.
Failure leaves accepted state, time and all cumulative budgets unchanged.

## Two different energy checks

The fluid-cell conservative equation gives, to roundoff under RK quadrature,

```
Delta sum(E_cell V_cell) + integrated outward boundary energy flux = 0.
```

This is local conservative accounting, not a whole-space accuracy certificate.
Independently compute

```
E_response = sum[(E_cell - B_cell,total²/(2mu0)) V_cell] + pi u^T K u/mu0
W_source = integral[-(2pi/mu0) (K u)_plasma · dpsi_external/dt] dt
response_residual = Delta E_response - W_source
```

The second budget contains actual gas energy plus induced self-field energy to
infinity. It excludes external self/mutual energy and is not coil-circuit input.
It is **not forced exact**: cell B energy, nodal magnetic energy, boundary traces,
CT fluxes and time quadrature differ at finite resolution. Their error must be
measured and decrease. Reports normalize it both to initial response energy and
to actual applied source work; a large initial gas energy must not hide a source
error. No energy correction is made to improve this diagnostic.

Physical assumptions remain those of a fixed fluid-filled ideal-gas chamber with
prescribed scalar eta. Moving plasma/vacuum interfaces, gas outflow, transport,
ionization, neutrals and Hall/kinetic validity remain separate requirements.
The original Stage 2 formation gate stays closed. See the measured report in
[CONSERVATIVE_VERIFICATION.md](CONSERVATIVE_VERIFICATION.md) for the historical
global method and [LOCAL_LIMITER_VERIFICATION.md](LOCAL_LIMITER_VERIFICATION.md)
for the local comparison and its remaining limitations.
