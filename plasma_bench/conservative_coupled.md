# Unsplit conservative MHD coupled to magnetic vacuum

`ConservativeVacuumMHD` evolves actual fluid-cell mass, momentum and total MHD
energy together with shared upwind constrained induction in SSP RK2. The
induced field outside the gas is completed by the existing exact discrete
elliptic/multipole solve at each stage. Prescribed external flux increments are
subtracted exactly when constructing induced state; analytic source power is
integrated with the same RK quadrature. No fluid exists outside the chamber.

`ConservativeState` energy includes cell magnetic energy, unlike `CoupledState`
from the historical split reference. States must be explicitly converted.
`ConservativeParameters` exposes gamma, physical eta, CFL, source time cap and
magnetic boundary controls. Reflecting gas ends/wall remain a verification
limitation; this is not a moving free-plasma interface or a formation claim.
Vacuum-side tangential fields provide the interface Maxwell stress and Poynting
flux. Shared normal CT fields and reflected normal gas velocities enforce the
fixed slip surface; interface work is not discarded by reflecting interior B.
Tangential field jumps use an explicit LLF interface electric dissipation in
both flux evolution and all adjacent Poynting transport; physical eta diagnostics
exclude this numerical effect and interface accuracy requires refinement.
The numerical interface eta uses D=a h/4 and the positive stiffness current,
not a one-sided derivative jump; its time-step restriction and independently
diagnosed positive numerical heat are explicit. `interface_diffusion=False`
retains the centered interface for comparison, including its boundary failures.
`spatial_order=1` exposes the piecewise-constant conservative/CT reference for
full-update admissibility studies; order 2 remains the default.
If a complete second-order Euler proposal is inadmissible, the optional default
`positivity_flux_limiting` selects conservative fluxes, CT electric fields,
metric sources and boundary-energy accounting from an admissible first-order
proposal. The default `flux_limiter="local"` replaces only shared entities
incident to troubled cells, rechecking neighbors until the mixed update is
admissible; `flux_limiting.md` documents its conditional termination argument.
The explicit `"global"` comparison retains the previous convex bisection.
Both SSP Euler stages are checked before their final convex average.
No cell is clamped or energetically repaired. A nonphysical first-order proposal
still raises atomically. Diagnostics record affected/troubled cell-stage counts,
maximum affected fraction, minimum high-order weight and iteration counts.
Per-cell cumulative masks permit independent locality checks and plots.
Counters commit only after the entire step validates. Physical/numerical eta and
source power at a given state are order-independent; outward energy is computed
from the actual mixed face fluxes. Limiting is order reduction and its activity
must be reported, not hidden as full second order.
Local conservation does not force the separately discretized full-space
response budget exact or guarantee its monotonic error under refinement.

Physical corner eta J uses induced-field stiffness and plasma-node dual areas;
external analytic vacuum fields have no volume current in the fluid. Its
Poynting flux enters total energy; physical Ohmic heating is diagnosed independently
and is not added twice. Rusanov/UCT diffusion is declared numerical dissipation,
not an added material resistivity or a separate tunable shock sensor.

Diagnostics distinguish exact fluid-cell energy-flux accounting from independently
computed gas-plus-full-space-induced-field response energy. The latter need not
cancel at finite grid/time resolution because finite-volume cell B energy and
variational magnetic energy differ. Its error must converge in source/boundary
studies; there is no global or local energy reset to force agreement. Source
work excludes external-field self/mutual or coil-circuit energy, as before.
Source cases report energy error relative to applied work as well as initial
response energy, so a large resting-gas energy cannot hide injection error.

All RK stages must remain physical before state, time and budgets commit. No
pressure/density floors, energy repairs or silent step fallback are used. The
historical split solver remains available for independent work-identity and
shock comparison. Original experimental Stage 2 remains gated.
