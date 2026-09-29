# External-engine result audit

Reopens every completed case under `outputs/engine_comparison/runs`, verifies
archive and driver/reference hashes, and recalculates continuum errors directly
from saved fields. Hall particle number and physical ion plus transverse magnetic
energy are reconstructed from individual macro-particle weights and momenta.
Vacuum energy includes electric and magnetic fields on their own staggered grids.
No energy residual is redistributed and no empirical trajectory acts as truth.
The Hall replay checks that the excluded axial background field stayed constant.
It also checks the initial electric field required by the continuum wave and
the particle integrator's first velocity half-step.
Every native energy-history step/time is checked, including maximum intermediate
energy variation. Hall native-history subtraction has an explicit text-precision
bound because background magnetic energy is much larger than the perturbation.
Independent endpoint particle/field archives avoid that cancellation. Native
input, diagnostic and measurement files are hashed in the final report.

Requires exactly the launcher's 17 cases, including three meshes (32/64/128)
for each Hall helicity and the vacuum wave, and native completed-time checks.
Time, loading, amplitude and floor controls are audited but excluded from the
mesh-refinement plot. The figure reports different models in separate panels,
and is explicitly labeled as verification rather than FRC formation.
`report.json` retains individual measurements and observed refinement orders.
Fixed-mesh differences between baseline, half and quarter particle time steps
are recorded separately from continuum errors to expose a spatial error floor.
Physical conclusions and migration recommendation belong in ENGINE_COMPARISON.
