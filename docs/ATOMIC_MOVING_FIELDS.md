# Moving variable-species electrostatic research component

## Limited numerical stress improves compression and current accuracy (2026-10-16)

Exact integration of J/E differences on the union of the saved physical meshes
confirms the original256→512 current difference:15.1346%. About88.20% of its
squared error is in |x|<9mm. The dominant discrepancy is an interior compressive
current structure, not the ion-support edge. The previous uniform sampling grid
was adequate for the spatial conclusion but missed narrow displaced cell jumps
in the much smaller same-grid time differences; corrected norms are below.

`MaterialSpecies(..., viscosity_limiter='minmod')` now removes resolved affine
compression from the numerical viscous jump using neighboring physical velocity
gradients. Q remains positive only in compression, and the identical stress
enters momentum and thermal work. No evolved density, current, energy or support
is clipped or repaired. Default `none` preserves the previous scheme; this is
an opt-in research improvement, not formation adoption. See both companions
for exact equations and the distinction from physical transport.

All original neutral/charged/plasma cases retain their masses, initial states,
vacuum and full durations. Complete actual20ns plasma comparisons are:

| Adjacent cells | Original J difference | Limited J difference | Original E difference | Limited E difference |
|---|---:|---:|---:|---:|
|128→256|20.5808%|21.3847%|1.56597%|1.47656%|
|256→512|15.1346%|13.0466%|0.912529%|0.705685%|
|512→1024|10.9910%|8.07805%|0.496145%|0.330304%|

Thus the limiter does not improve every coarse comparison. At1024 cells its
electron/ion fronts are±54.065224/±10.043085mm, numerical heating7.70463e-9J/m²,
and max relative total-energy error2.07849e-10. Density electron L1 difference
512→1024 is0.423637% of the initial number. All actual thermal energies remain
positive. The neutral front error is still1.19915mm; its slow convergence is
unresolved. Current peaks remain spatially sensitive despite improved fields.

Same-grid half-step J differences with exact union integration are0.174862%,
0.0545214%,0.0196930% for the limited32/64/128 cases. For the already archived
original-stress half-step cases they are0.175964%,0.0242791%,0.0287379%.
These supersede the earlier uniform-sampling J time norms, which missed tiny
displaced density discontinuities. E time differences stay below4e-6 relative.
The time differences are still much smaller than spatial differences.

Two independent physical baselines evaluate the actual change:

- Smooth neutral quadratic-pressure compression/rebound: scalar energy-integral
  quadrature sets the full14.3312877µs duration and minimum scale0.5. Original
  stress produces endpoint position errors1.03627/.465801/.218267mm at32/64/128;
  limited errors8.84e-12/6.20e-13/4.03e-14m. Its numerical heating is~1e-16J/m²
  instead of~9.08e-4J/m² at128 cells. Global conservation alone did not expose
  the original spurious heating and trajectory error.
- Existing analytic gamma5/3 Sod solution, full3.860718µs, finite free neutral
  slab and exact exterior vacuum. Central normalized density L1 errors for
  limited64/128/256 decrease.0119667/.00665079/.00353818 (original finest.00392436).
  Finest velocity/pressure L1 decrease.0102663/.00528047→.00803898/.00478187.
  Small rarefaction-head overshoot persists; it is not claimed monotone.
  Finest max relative energy error1.04647e-7, momentum change6.05e-21kg/(m s).

All runs below were collected EXIT0; their actual histories and frozen sources
are retained. No unit suite or generic verification framework was introduced.

| Output under outputs/ | Handle | Seconds | Source prefix |
|---|---:|---:|---|
|species_material_limited_support|23277|4.515|aca2b117|
|species_material_limited_refined|78703|24.982|aca2b117|
|species_material_limited_n1024|6725|78.138|aca2b117|
|species_material_original_n1024|98820|75.397|aca2b117|
|species_material_limited_halfstep|73375|9.089|aca2b117|
|species_material_smooth_compression|46438|3.384|824e580f|
|species_material_shock|18095|1.935|dfdc0f8c|

`outputs/species_material_current_localization/analysis.py` reproduces the exact
saved-field comparisons and their PNG/PDF, including corrected historical time
norms. Compression, shock and current figures were visually inspected. Bead
frc-transport.5.5.1 and the full goal remain OPEN. Next resolve the narrow
interior current structures and low-mass fronts with improved spatial treatment,
then couple physical changing populations and their momentum/energy across
distinct species supports. This remains a reaction-disabled planar component;
fixed numbers are not a substitute for atomic reactions or full formation.

## Separate species supports remove the refined thermal failure (2026-10-16)

New `plasma_bench/species_material` advances each species on its own physical
material support, with consistent nodal mass, three velocities and cell internal
energy. The exact represented Gauss field uses the union of the species cells
and the empty exterior. Electric work enters actual nodal kinetic energy rather
than an ALE numerical-current thermal source. Read the module/example companions
for the complete equations, boundaries and disclosed numerical viscosity.
The full variable-population/reaction and axisymmetric goals are unchanged:
this reaction-disabled transport component is NOT adopted by the formation solver.

All original neutral/charged/plasma physical cases and durations were retained.
Physical masses, original32/64/128 grids, then128/256/512 grids all complete;
the electron–ion request is20ns on every grid. No clipping, density floor,
exterior matter, mass deletion, shortened duration or thermal-energy reset.
Every accepted state and exact source is archived in these completed EXIT0 runs:

| Output under outputs/ | Handle | Seconds | Source prefix |
|---|---:|---:|---|
|species_material_support|20816|4.326|6ce08375|
|species_material_cosine_support|17772|4.806|e72eaed6|
|species_material_cosine_viscous_support|71020|4.407|18e7b306|
|species_material_cosine_viscous_halfstep|48102|8.337|18e7b306|
|species_material_cosine_viscous_refined|88244|25.653|35bd273f|

Uniform material particles resolve the low-mass rarefaction front slowly.
Cosine-spaced initial labels increase surface resolution without changing
the physical initial density or inventory. The inviscid version has substantial
compressive oscillations. Explicit numerical compression stress C1=.5,C2=2
damps those oscillations with matched mechanical work and ion/electron thermal
heating; it is numerical viscosity, not a transport coefficient inferred for a lab.
No-viscosity attempts remain archived. Charged affine expansion has no numerical
viscous heating and reproduces the independent energy-integral reference to
time/roundoff accuracy; its512-cell front error is8.95e-16m.

For cosine labels with the declared numerical stress:

| Cells | Neutral number L1 / initial number | Neutral front error | Plasma max relative energy error |
|---|---:|---:|---:|
|32|0.0260762|6.23687mm|2.70232e-8|
|64|0.0141323|4.57148mm|6.74310e-9|
|128|0.00709053|3.30380mm|1.77860e-9|
|256|0.00347507|2.36704mm|1.09294e-9|
|512|0.00172143|1.68713mm|4.65887e-10|

All plasma histories conserve particle amounts exactly by their material
representation and keep momentum near roundoff. At512 cells the electron fronts
are±53.294637mm and ion fronts±10.038820mm. Adjacent128→256 and256→512 electric
field differences at20ns fall1.56591%→0.912610%; current differences fall
20.5750%→15.1339%, and electron-density L1 fractions1.4680%→0.879826%.
Thus the small current-carrying structure and leading front still need spatial
resolution; these results do not justify a converged physical edge claim.
The measured exterior field energy is retained even where every species is empty.

Halving time limits at fixed32/64/128 meshes changes finalJ by at most
4.12770e-5 relative (128 cells4.39032e-6) and E by at most3.49246e-6.
These time differences are far smaller than the spatial differences. Exact
same-grid endpoint comparisons are saved in the half-step output. Fine plots
were inspected and their PDFs saved. No unit tests or new generic verification
framework were used; this is actual replacement transport/front dynamics.

Beadfrc-transport.5.5.1 remains IN_PROGRESS. Next improve the low-mass front and
local current spatial accuracy, then conservatively couple changing atomic
populations, reaction momentum/energy and collisions across distinct species
supports. Do not silently keep fixed species numbers when making that adoption.
Transverse Maxwell fields, cylindrical geometry and physical wall/end/sheath
conditions remain full-goal requirements. Earlier ALE failures below are preserved.

This component advances actual species numbers, three momenta and total energies
on moving axial control volumes. It includes relative species fluxes, a Gauss
field from actual charge, and electric work compatible with the moving-volume
continuity equation. Field energy is integrated over the material region AND
the exactly empty exterior out to fixed remote field boundaries. It neither
enforces neutrality nor resets carrier numbers, temperature, current or energy.

It is **not ready for formation adoption**. The outer-front trajectory remains
inaccurate, and refined electron–ion runs reject negative thermal energy in tiny
numerically diffused ion tails. Exact emptiness outside a numerical front does
not prove that the front is physically correct. The existing axisymmetric
formation solver and its live source closure have not changed.

Equations, assumptions and numerical diffusion are in the atomic_moving,
species_ale_flux and species_reconstruction companions. Independent neutral
rarefaction and warm charged affine-expansion references are derived in
`examples/atomic_moving.md`. The latter uses an independent scalar energy-integral
quadrature, not another integration of the discrete solver.

## Preserved complete attempts

All four runs are terminal EXIT1 because at least one requested case failed.
Each preserves every accepted state, source archive and the attempted step.

| Output directory under outputs/ | Process | Seconds | Source SHA256 |
|---|---:|---:|---|
| atomic_moving | 47543 | 2.085 | bf237432abe5d26da929c816f5ba2aa80dcfc6295dd681c8ce97886bef78e650 |
| atomic_moving_cell_front | 61410 | 19.327 | 9f7aa13490ec3782e922f751ca1eae8a32850684cd0d023222d6c0eb3cd3576c |
| atomic_moving_affine | 46366 | 9.279 | e440ab212edd388d18f931f7937550ad653b8f87fe02212e476d4c789dd3f7a3 |
| atomic_moving_conservative | 83790 | 53.439 | c96f8697718d8fff478e95f2aca00c635c0655608c96f9247ea8f823ed354e7e |

The first attempt incorrectly tried to determine a vacuum-front sound speed
from an exactly empty reconstructed trace. Posing the vacuum Riemann problem
from the adjacent occupied cell removed that undefined operation. Barycentric
mesh motion then stretched the outermost cell across the charge-separation
region: the finest failed electron cell spanned approximately -19.84 to -9.84 mm.
Its reconstructed rejected stage is preserved in atomic_moving_cell_front.

Affine coordinate motion now keeps relative cell widths fixed while physical
species velocities remain independent. This alone did not fix thermal failures.
Conservative reconstruction now preserves the mean number/momentum/total energy;
a common admissibility slope factor acts only on face traces. This replaced
independent primitive reconstruction, which could misrepresent kinetic energy.
No evolved state was clipped, floored or repaired. This is numerical limiting,
not a claim of a positivity theorem for the complete moving electric system.

## Latest measured results

All neutral and charged-reference cases now finish. Neutral duration is
9.850823 us; the charged reference reaches scale s=1.5 at 4.446830 us.

| Cells | Neutral number L1 / initial number | Neutral front error | Charged number L1 / initial number | Charged front error |
|---|---:|---:|---:|---:|
| 32 | 0.0367999 | 1.54072 mm | 0.0624625 | 19.9082 mm |
| 64 | 0.0208654 | 1.29066 mm | 0.0494690 | 16.1442 mm |
| 128 | 0.0140413 | 5.25630 mm | 0.0415637 | 13.3992 mm |

Density improves, but the support/front does not satisfy convergence acceptance.
In the charged case the true final half-width is 37.5 mm; the finest numerical
front is 50.8992 mm. Finest field error relative to the exterior field is 1.424%.
Maximum relative total-energy error falls 5.594e-9,1.044e-9,2.171e-10. These small
energy residuals do not excuse the large front error. Particle inventory and
total momentum stay near roundoff; the field includes exterior energy throughout.

The neutralized electron–ion case starts with a 20 mm-wide slab, n_i=n_e=1e14/m³,
Ti/Te=0.2/2 eV and zero E. The requested duration remains 20 ns on every grid:

| Cells | Accepted steps | Last accepted time | Outcome |
|---|---:|---:|---|
| 32 | 290 | 20 ns | Complete; front ±69.9448 mm, peak E443.792 V/m |
| 64 | 391 | 13.4072 ns | Rejected step392, negative ion thermal energy |
| 128 | 252 | 3.65695 ns | Rejected step253, negative ion thermal energy |

The reconstructed actual failed RK4 stages are saved in `rejected_stage_n64.npz`,
`rejected_stage_n128.npz` and `rejected_stages.json`. Both are stage4, cell0, ion:
number/Nref=1.71537e-27 and1.12953e-24; thermal energy/Nref=-6.25320e-38 and
-1.33521e-34 eV. These are numerical ion tails within the outer electron support,
not a reason to delete particles or allow negative temperatures. The last
accepted charge and energy budgets remain small, illustrating the local failure.

The original final-state plot did not label different stopping times. It remains
archived; `moving_species_annotated.png/.pdf` explicitly labels actual endpoints
and completion status. It was generated from the same saved arrays without
rerunning physics. The current example changes those future plot labels only.
Do not interpret the three plasma endpoint curves as an equal-time comparison.

## Concrete next work

Replace repeated cell-mean vacuum-fan speed estimates with a consistent evolved
front state/pressure treatment, and distinguish individual species support when
electrons outrun ions. Resolve local thermal admissibility through the numerical
equations rather than a density cutoff or energy reset. Preserve the same three
physical cases, requested times, physical masses and independent references.
Require converging front positions, complete refined electron–ion histories,
equal-time field/current comparisons and the independent global budgets before
coupling this component to the formation geometry. Spatial atomic reactions,
transverse fields, radial geometry and real boundaries remain original-goal work.
