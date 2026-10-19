# Unsplit conservative MHD coupled to magnetic vacuum

`ConservativeVacuumMHD` evolves actual fluid-cell mass, momentum and total MHD
energy together with shared upwind constrained induction in SSP RK2 or RK3. The
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
Every SSP Euler proposal is checked before its convex combination.
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

Experimental `energy_quadrature="staggered"` uses positive face-square energy
consistently in construction, pressure, admissibility and budgets.
`boundary_poynting="variational"` uses the natural discrete vacuum traction
and the actual selected CT electric field for boundary work. Both are explicit
controls during verification; the historical cell/Riemann choices remain
selectable. No evolved state is converted to change its thermal energy.

`time_order=3` selects the nonautonomous three-stage SSP method with RHS times
`t, t+dt, t+dt/2` and budget weights `1/6, 1/6, 2/3`. Stage completion subtracts
the external field at that stage's actual time. The historical RK2 remains
the default comparison. All stages and their limiter counters commit atomically;
`euler_stages` records the accepted count instead of assuming two per step.
Neither method preserves a general quadratic magnetic-energy invariant exactly.
`advance(..., on_step=callback)` optionally observes each accepted step for
long-run progress; the callback must not mutate the solver. It is not invoked
for rejected proposals and does not alter timestep selection.

## Opt-in dual energy and timestep halving (2026-10-19)
`dual_energy=True` adds an advected entropy density ρA, with A = p/ρ^γ (see
`dual_energy.md`). It addresses the low-β positivity limit: in cells where the
gas energy is a tiny fraction of the magnetic energy, pressure recovered from
total energy carries a relative error of roughly 1/β.

**Defaults:**
- `dual_energy_mode="pressure"`: the conserved total energy is never changed.
- `dual_energy_threshold=0.1`.
- `dual_energy_shock_threshold=0.2`.

**How it works:**
- Cold cells take their pressure from ρA.
- The limiter exempts cold cells from the internal-energy check only.
- Warm cells resynchronize ρA from total energy after every step.
- ρA advances by one first-order forward-Euler step per time step. Its physical
  Ohmic heating source is included; the numerical interface heat is not.
- `reset` mode instead sets the total energy of cold cells to the entropy
  pressure and counts the injected energy (`dual_energy_change_j`).
- In pressure mode, `dual_energy_cold_discrepancy_j` reports the current
  (E − K − M) minus the entropy thermal energy summed over cold cells.

`max_step_halvings=N` (default 0) retries an inadmissible step from the
accepted state with dt/2, up to N times, and counts the retries in
`step_retries`. It is timestep control, not a state repair; every stage still
validates atomically.

Both options default off. With them off, results are bit-identical to the
previous code on the standing-shock and coil-pulse fixtures (2026-10-19
comparison against commit ee9621b).

`local_trace_bounds=True` passes the reconstruction's neighbour-minimum bounds
for density and thermal pressure (see `mhd_reconstruction.md`). Off by default.

`minimum_density` (kg/m³; default None) enables the declared gap closure in
`density_floor.md`, with `floor_temperature_ev` and `floor_ion_mass`. It runs
after the MHD step and before commit, inside the same atomic step attempt.

`resistivity_model="spitzer"` replaces the uniform `resistivity` with the NRL
perpendicular Spitzer value, 1.03e-4 lnΛ T_e^-1.5 Ω·m.
- **Inputs:** T_e = p/(2 n e), with n = ρ/`spitzer_ion_mass` (Z = 1, T_e = T_i),
  lnΛ = `spitzer_coulomb_log` (default 10), and T_e floored at
  `spitzer_min_te_ev` for evaluating η only.
- **Placement:** cell values are averaged to nodes with the coupling mesh's
  cell→node weights.
- **Consistency:** the same node η enters the electric field, the Ohmic power,
  the dual-energy heating source and the diffusion time-step limit, so the
  energy accounting stays exact.
- **Plumbing check:** with an ion mass chosen so the Spitzer value equals the
  coil-pulse fixture's uniform η at its initial temperature, states agree to
  2e-6 at t = 0.01 and then drift only as heating changes T. The fluid energy
  residual stays at roundoff.
