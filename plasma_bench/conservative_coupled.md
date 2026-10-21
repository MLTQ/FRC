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

`chodura_coefficient=c > 0` adds the Chodura (1975) anomalous resistivity used
by theta-pinch formation codes: η_ch = c m_e ω_pi/(n e²) [1 − exp(−v_d/(f c_s))].
- **Inputs:** n and c_s are averaged from cells to nodes; v_d = |J|/(n e) at
  the node; ω_pi uses `spitzer_ion_mass`; `chodura_f` defaults to 3.
- **Where it acts:** where the electron drift exceeds about f c_s, i.e. in thin
  current sheets such as the end reconnection sites. It grows as n^-1/2 in
  tenuous plasma.
- **Status:** a declared model, not a fit, with c ≈ 0.1 from the literature.
- **Check:** recomputing it by hand at the node of largest current reproduces
  the solver's value exactly, and energy closes to roundoff.

`interface_speed_cap` (m/s; default None) passes that bound to the interface
numerical resistivity (see `vacuum_trace.md`).

`open_ends=True` makes the gas ends outflow-only (see `vacuum_flux.md`). The
mass leaving is tracked as `outflow_mass_kg` and included in
`mass_relative_change`. It is computed exactly from the update, since interior
fluxes telescope.

`open_line_flux_limit=f > 0` enables the declared conduction loss in
`open_line_loss.md`, applied after the closure step and counted in both
residuals.

`dual_energy_criterion="entropy"` (default `"total"`) selects the cold-cell test
described in `dual_energy.md`. Use it with `dual_energy_mode="reset"`: the
pressure-mode/total-criterion pair can turn a hot, weakly magnetized cell into
an unbounded energy reservoir.

Reset-mode corrections (`dual_energy_change_j`) are now a declared source in
both the fluid and the response energy residuals. In pressure mode they are
zero. Before this change, the reset-mode fluid residual equalled the reset sum.

`analytic_source_traces=True` builds the wall and end tangential traces from
the differenced *induced* flux plus the analytic coil field (`vacuum_trace.md`).
A coil filament inside the one-sided trace stencil then cannot corrupt them.

`magnetic_energy_flux` (default `"face"`, the historical scheme) selects how
magnetic energy is transported in the total-energy flux (`vacuum_flux.md`):
- `"dissipative"`: the LLF magnetic dissipation is carried as the Poynting flux
  of the UCT dissipative field.
- `"ct"`: all magnetic energy is carried as the Poynting flux of the full CT
  electric field, which is the field that changes ψ.
- `ct_dissipation_heat_j` estimates ΣE_diss·J from the high-order fluxes.

**Measured in NUCTE-III case 6** on the compatible path, by single-step
decomposition of the drift between cold-cell total and entropy energies (the
drift reset mode corrects), at 2.5 and 3.5 µs:

| Mode | Drift magnitude at 2.5 µs (W) | Drift magnitude at 3.5 µs (W) | Wall-row drift |
|---|---|---|---|
| `face` | 8.3e9 | 1.1e10 | 4.2e9 |
| `dissipative` | 4.4e9 | 4.0e9 | — |
| `ct` | 2.9e9 | 1.4e9 | 1.4e8 (30× smaller) |

On the historical cell/Riemann path, `dissipative` made no difference to the
reservoir failure. That failure was a criterion error, not a flux error.

`gap_resistivity` (Ω·m; default 0) and `gap_density` (kg/m³) add a declared
resistivity η_gap/(1 + (ρ/ρ_gap)⁴) at nodes. ρ is averaged from cells, so the
term acts on near-floor cells: the wall–piston gap and the un-pre-ionized gas
beyond the coil.

**Why:**
- The declared density floor is a fully ionized 2 eV plasma. Spitzer resistivity
  does not depend on density, so this tenuous floor plasma conducts as well as
  the real plasma and freezes the bias flux.
- Beyond the pre-ionized zone, that frozen flux keeps the reversed field reaching
  the tube ends. In NUCTE-III case 1 and in a uniform 12" coil with D2 this
  blocks end reconnection.
- Un-ionized gas is an insulator, and the gap is near-vacuum. A large η makes
  both magnetically transparent.

**Accounting:** the same node η enters the electric field, the Ohmic power and
entropy heating, and the diffusion time-step limit, so energy accounting is
unchanged.

**Cost:**
- Diffusion across a length L takes μ0L²/η_gap.
- Making 15 cm transparent within a 3 µs rise needs η_gap ≳ 0.05 Ω·m, which
  shrinks the explicit time step about 5×.
- 1e-2 Ω·m clears a few-cm gap in about 0.3 µs.

## Annular chamber (coaxial inner tube)
An `MHDMesh` with `inner_radius > 0` gives a gas annulus between an insulating
inner tube and the outer wall. The inside of the tube is field-only vacuum on the
coupling grid, which still starts at the axis.
- Source filaments may lie inside the tube (r < inner_radius) or outside the
  wall, but not in the gas.
- The induced state is stored at field rows `r_offset..`.
- The inner wall reflects gas and is magnetically transparent.
- Traces, Poynting flux and quadrature are described in vacuum_trace.md,
  vacuum_energy.md, vacuum_flux.md and magnetic_quadrature.md.

**Verification** (`examples/annulus_verification.py`):
- **Static equilibrium:** a uniform plasma in a uniform field stays at rest,
  with max speed 6e-10 m/s against v_A 9e4 m/s and density change 6e-14.
- **Resistive flux penetration:** flux enters the inner vacuum only through a
  resistive annulus (r 4–10 cm, η = 2e-3 Ω m, applied field 0.05 → 0.25 T in
  2 µs). It is compared at the midplane with an independent 1D solution in
  x = r²/2.

  | nr (annulus) | max ψ error / max induced ψ (1, 2, 3, 4 µs) | inner-tube Bz error at 1 µs | energy residual |
  |---|---|---|---|
  | 12 | 2.9e-3, 7.7e-4, 1.5e-3, 1.9e-3 | 1.4e-3 | 2e-8 |
  | 24 | 7.1e-4, 1.9e-4, 3.8e-4, 5.0e-4 | 3.5e-4 | 6e-10 |
  | 48 | 1.8e-4, 4.9e-5, 9.8e-5, 1.5e-4 | 8.8e-5 | 2e-11 |

  The error falls by 4.0 per doubling through 3 µs, which is second order. At
  4 µs it falls by 3.8 and then 3.2, where the induced flux is small and the
  residual motion of the heavy fixture plasma (0.03 m/s) begins to show. The 1D
  reference's own discretization error is 8e-8.
- **Axis case unchanged:** `inner_radius = 0` reproduces the axis runs bit for
  bit (a NUCTE formation snapshot checked against the previous commit).

## 2.5D: rotation, toroidal field and the Hall term (2026-10-20)
`toroidal=True` adds two conserved cell variables: L = ρ⟨r²⟩Ω (the average of
r ρ v_θ) and g = B_θ/r (weighted by r dr). Both start from the state the caller
gives (`mhd_state.with_toroidal`). Default off: the poloidal-only path is
bit-identical to the previous code. This was checked on the NUCTE default and on
frc12_uniform_h2 to 2 µs: all 20 snapshots, velocities and history rows equal,
and the 2 µs frc12 states equal `outputs/fine_uniform_h2`. With `toroidal=True`
and no seed or Hall term, v_θ = B_θ = 0 stay exactly zero and the poloidal
state is unchanged bit for bit.

**Equations.** Single-fluid MHD in r–z with all three components:
- angular momentum ∂(rρv_θ)/∂t + ∇·[r(ρ v_θ v − B_θ B/μ0)] = 0;
- Faraday's law for B_θ on each r–z cell (`toroidal_ohm.md`), with
  E = −v×B + ηJ (+ J×B/(ne) with Hall);
- the radial momentum source gains ρv_θ²/r − B_θ²/(μ0 r), and the total
  pressure and energy include B_θ²/(2μ0) and ρv_θ²/2.

The Rusanov fluxes use `toroidal_flux.physical_flux`, which is regular on the
axis, and the fast speeds include B_θ (UCT, time step, interface resistivity).

**Boundaries.** The glass wall, the ends and an inner tube are insulating, and
the vacuum outside is exact. v_θ slips.
- **Default.** No current links the vacuum, so B_θ = 0 on every gas boundary
  and no B_θ energy crosses one.
- **`toroidal_boundary=rbt`** (annular chambers only) is a callable
  t → r B_θ [T m] of the vacuum beyond every gas boundary, μ0 I(t)/(2π) for a
  current I threading the vessel (a rod on the axis inside the inner tube,
  frc-ibl.9.1).
  - The boundary value enters the reflections, the currents and the boundary
    Poynting flux.
  - The B_θ energy that this flux delivers is counted separately as
    `toroidal_boundary_work_j`, and as a source in the response residual.
  - Verified in `examples/toroidal_verification.md`, case `rod`. B_θ diffuses
    into a resistive annulus from both walls. The midplane rB_θ matches an
    independent 1D solution to 4.6e-3, 1.2e-3 and 3.2e-4 of its final value at
    12, 24 and 48 cells. The boundary work equals ΔW_θ + Ohmic heat to
    1.4e-5, 1.2e-6 and 1.3e-7.

**Resistivity.** η J_z and η J_r use the same node η (Spitzer, Chodura, gap).
The Chodura drift uses |J| including J_r and J_z. The paired heat
Σ η J² V enters `ohmic_heat` and the dual-energy entropy source.

**Hall term (`hall=True`, DECLARED).** E_H = J×B/(n e), with n = ρ/
`hall_ion_mass` (default `spitzer_ion_mass`). It is built skew-adjoint, so the
discrete Hall field does no work. `hall_work_j` is the time integral of
Σ J·E_H, reported (measured ≤ 1e-16 J). The electron-pressure term ∇p_e/(n e),
which would add Biermann-battery B_θ, is omitted; so are electron inertia and
two-temperature effects.

Requirements and limits:
- **RK3 only.** The term needs `time_order=3`: whistlers are undamped, and
  only RK3's stability region contains part of the imaginary axis.
- **Time step.** The CFL sum adds cfl ω/(f √3), with
  ω = 4 max(c|B|)/μ0 (1/dr² + 1/dz²), c = w/(n e) and f =
  `hall_whistler_fraction` (default 0.5). In the whistler test, f = 1.5 is
  stable and f = 2 grows, so the ω bound is about 1.7× conservative and the
  default keeps about 3.5× margin. `hall_limited_steps` counts the steps where
  this term dominates.
- **Density cutoff (DECLARED).** `hall_cutoff_density` ρ_c applies the weight
  w = 1/(1+(ρ_c/ρ)⁴), so c ≤ 0.57/(n_c e). Near-vacuum floor cells (1% of fill)
  have d_i 10× the fill value. Without the cutoff they would cut the whistler
  step about 100×, and their whistler physics is not meaningful.
  `hall_cutoff_cell_steps` sums the cells with w < 0.5 over steps, reported
  separately like the gap resistivity. The formation driver sets ρ_c to 0.1 of
  the fill density.

**Verification.** `examples/toroidal_verification.md` compares against exact
solutions, using the formation settings (staggered energy, variational boundary
Poynting, RK3, ct, local trace bounds):

| Case | Measure | 8 | 16 | 32 | 64 |
|---|---|---|---|---|---|
| Torsional Alfvén standing wave, 1 period | B_θ L1/amplitude | 3.0e-3 | 8.3e-4 | 2.0e-4 | 5.1e-5 |
| | frequency error | 3.6e-3 | 1.6e-3 | 4.5e-4 | 1.0e-4 |
| | toroidal energy change | −14% | +1.0% | +0.38% | +0.062% |
| Isothermal rigid rotor, 3 sound crossings | density L1 | 7.3e-4 | 1.1e-4 | 1.6e-5 | 2.2e-6 |
| | r ρ v_θ L1 | 2.3e-4 | 3.0e-5 | 3.9e-6 | 6e-7 |
| Alfvénic toroidal state (v = B/√μ0ρ), Mach ~1, 2 crossings | B_θ/r L1 | 1.7e-2 | 2.6e-3 | 5.3e-4 | 1.0e-4 |
| | pressure L1 | 1.1e-2 | 1.7e-3 | 3.1e-4 | 6.2e-5 |

The toroidal energy of the torsional wave rises by the listed amount in the
first period, then decays (at 16 cells: −1.7% by 6 periods). The offset
converges and is the same in the `face` and `ct` modes.

Whistler (EMHD packet, immobile ions, exact Hankel–Fourier propagator; cells per
wavelength λ/dz ≈ 7, 14, 29):

| Whistler measure | 36×94 | 72×187 | 144×374 |
|---|---|---|---|
| frequency error (fitted scale on D) | −6.8% | −1.75% | −0.44% |
| B_θ relative L2 (energy norm), 1 period | 48% | 12.5% | 3.2% |
| same after the frequency fit | 15% | 3.7% | 0.99% |
| induced poloidal flux relative L2, T/4 | 4.7% | 1.3% | 0.34% |
| share of the B_θ error within w_r/4 of the axis | 4.3% | 1.6% | 2.0% |

The frequency error matches the −(k dz)²/12 dispersion of the centred
operators: −6.3% at λ/7.

Fluid energy residuals are at roundoff in every case. In the whistler case the
response-energy residual (1e-13 J) is below 1e-4 of the 7.5e-10 J exchanged
between B_θ and the poloidal field.

**First formation run with Hall** (frc12_uniform_h2, cutoff 0.1 of the
fill): see `docs/NUCTE_BENCHMARK.md` gap 5.
- **Hall limits the step rarely.** 349 of 181,930 steps. The extra cost comes
  from the additive whistler term: 1.68× the steps and 1.31× per step.
- **Cutoff.** On average 1,561 of 3,456 cells sat below the cutoff weight 0.5,
  mostly the un-pre-ionized ends.
- **Hall work.** −3e-17 J.
- **Energy.** The response residual is 9.6e-5, against 8.9e-5 in 2D.

**Angular momentum.** Without a threading current, B_θ = 0 on every
boundary, so no magnetic torque crosses one, and the slip boundaries exert
none. Total angular momentum is
therefore conserved to roundoff. In a z-symmetric formation the Hall term makes
B_θ odd in z and v_θ even in z, changing sign with r: it redistributes rotation
but cannot spin the plasma up as a whole. Net ion spin-up (end shorting,
particle loss) is not modelled.

**Known limitation.** Near the axis the skew Hall pairing has an O(h²/r) error,
50% of J_r at the first node row (`toroidal_ohm.md`). Field norms converge; the
B_θ/r norm is dominated by the axis cells.
