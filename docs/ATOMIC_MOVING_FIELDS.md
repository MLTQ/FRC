# Moving variable-species electrostatic research component

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
