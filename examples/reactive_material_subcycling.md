# Full charged argon with explicit reaction intervals

`--electron-cells` optionally selects a distinct initial electron density and
velocity mesh. Its default is `--cells`; all atomic species retain the latter
initial/newborn target. Each species still has its physical mass and complete
support, and return retains every inherited knot and density cut. The common
plasma-phase cap is `.04*32/max(cells,electron_cells)`, recorded explicitly and
unchanged by source enrichment. Compare with the existing equal-resolution
case at that same cap before claiming that cheaper heavy-species meshes suffice.
This targets the measured dominant interior electron-current error; it is not
a fluid closure change, a population cutoff or an accuracy guarantee.

`--electron-spacing uniform` changes only the initial electron grid from the
default cosine spacing. Every cell still represents the same constant initial
density/temperature on[-10,10]mm; amounts are density times actual cell width.
Heavy grids remain cosine-spaced. The initial electron count and phase cap are
unchanged. At128cells, centre spacing changes245.412 to156.25um, while edge
spacing changes3.012 to156.25um. This shifts resolution toward the measured
interior current discrepancy and relaxes the initial tiny-edge acoustic bound;
it may worsen the vacuum-fan approximation and must be compared over full20ns.
It removes no evolved knot, density cut, population or physical endpoint, and
changes no masses, equations or transport bound. Reports record actual spacing.
Earlier references described as uniform128 mean128cells for every species;
those runs use cosine spacing and are not this new uniform-spacing experiment.

Uniform spacing explicitly supplies every initial electron face as a velocity
knot. The first experiment omitted this argument, so the constructor selected
102cosine-target velocity cells on128uniform density cells. Its complete20ns
output is preserved as a mixed-mesh attempt, not a128-uniform-velocity result.
The corrected run starts with128actual uniform velocity cells and preserves
them through return. The summary now also records actual initial velocity-cell
counts. Saved initial indices, not only CLI counts, define the mesh.

Retains the full20ns neutral/ion/electron ground-ladder argon problem and its
physical masses,1e14/m³ initial component densities,10mm half-width and
0.1/0.2/2eV temperatures. The precise geometry and vector-assembled source are
unchanged. Physical-width coarsening is absent. `--coalesce-ulps` now defaults
to zero, retaining every distinct density cut. The earlier recorded10ns/5ns
attempts used8-ULP internal cleanup; that remains an explicit comparison option
and still protects every species endpoint. Existing atomic/closure limits remain.

Reaction return now preserves existing velocity knots and records the actual
per-block velocity-cell counts. The CLI cell count sets the initial mesh and
new-species target; inherited knots and new support endpoints can increase it.
Saved indices define the actual mesh. No knot is removed to enforce a quota.

`--source-step-ns` chooses a disclosed source interval (default10ns), not the
motion timestep. Within it the material equations/fields use actual RK4 steps
at half their current acoustic/plasma bound. Every accepted internal dt and
physical budget is retained separately. Reaction states can change rapidly
through the motion, so adequate source timing must be demonstrated by complete
10/5/2.5ns or finer comparisons; no convergence claim follows from completing
the first chosen interval. No implicit transport or enlarged physical mass.

Each completed coupled interval saves the full precise state and reservoir
budgets. Internal transport stages save times/budgets, not full geometry
snapshots; they are explicitly split-stage records, not completed coupled states.
Source/input identities, initial state and actual internal step sizes retain
reproducible evolution. A rejected source/transport interval retains its last
source and internal diagnostics but does not advance completed coupled time.
No clipping, hidden retry or shortened success. Disk reserve remains256MiB.
Rejected-source archives now include `velocity_knots_species*_decimal_m` so
the actual inherited projection can be reconstructed. The first retained-knot
5ns attempt predates this extra metadata; its complete source and transport
step records remain, but its failing transfer requires evolution replay.

Use `python -m examples.reactive_material_subcycling FRESH_OUTPUT
--source-step-ns 10 --cells 32`. This is an opt-in coupled research calculation,
not a replacement for the required spatial, source-time and transport-time
convergence, physical collision/field/wall integration or FRC formation evidence.

`--projection-method bounded_moments` selects the explicit conservative bounded
momentum/energy transfer. The default remains `galerkin`. Reports preserve the
method, both cumulative numerical thermal/redistribution diagnostics and the
last source return's local limiter activity. This is a comparison of numerical
transfer methods; it does not change the atomic source, masses or resolution.
NumPy scalar diagnostics in the nested projection record are serialized as
ordinary JSON numbers, preserving their values.

`--projection-method constrained_velocity` selects the volume-norm velocity
projection with physical block momentum and each fine cell total energy
constrained. Source and transport physics are unchanged; barrier/constraint
activity is a numerical approximation, not an adoption criterion.

Print accepted internal transport progress every100 steps. Terminal interval
stdout is compact; the complete per-block projection record stays in JSON.
No accepted states or actual timestep records are omitted from saved results.

`--extend-velocity` explicitly continues inherited edge basis functions across
new support material as documented in velocity_extension.md. Every old velocity
knot, density cell and physical support endpoint remains. This option and the
projection method are distinct; neither changes the source interval or declares
the local thermal/current errors converged. Saved states include extension flags.
# Newborn mesh comparison

The report records imported geometry arithmetic and bit precision. Select the
optional 192-bit backend with `FRC_GEOMETRY_BACKEND=mpfr` before process start;
install `plasma-bench[precision]`. Decimal remains default. Saved physical
arrays and full-duration acceptance are independent of the selected backend.

`--newborn-spacing-fraction` explicitly selects the target-relative initial
mesh rule in `velocity_selection.md`; default zero retains the earlier rule.
Reports record the fraction and every block's actual dimension. No existing
velocity knot or density cell is removed. The option does not change the full
20 ns target or turn prefix progress into a successful trajectory.

`--projection-method variance` evaluates the explicit conservative velocity-shape
and mixing-heat transfer described in `plasma_bench/variance_projection.md`.
The original full duration and physical controls remain. This option is
experimental; it must pass local accuracy and complete trajectory comparisons.

`--velocity-refinement-rounds` requests the explicit existing-face mechanical
enrichment in `velocity_refinement.md` at each source return, default zero.
`--maximum-refined-velocity-cells` limits additions by whole rounds (default
256), retaining every inherited node even above that capacity. Reports record
both controls and every actual per-block addition/capacity limit. Density cuts,
source quadrature, populations and full duration stay unchanged. Actual transport
steps obey the enriched acoustic bound. This changes spatial resolution and
must be compared as such, not presented as a pure source-time refinement.

`--refine-heavy-species` restricts additions to the network's heavy atomic
species and records their actual indices. Electrons still evolve on their
inherited space with physical mass/current; no carrier is omitted. This isolates
the measured heavy-ion heating problem. The all-species enriched run increased
cost and changed current by 14.67% but barely improved its discrepancy from the
existing fine numerical reference, so late electron enrichment is not treated
as a current-accuracy fix. Default selection remains all species when refinement
is explicitly enabled; default refinement itself remains off.

`--normalization-shift-bits` selects the explicit binary coordinate scaling
of reaction return. Default0 is unchanged;128 raises stored relative N/U by
2**128 and lowers their common physical scale by the same amount. All physical
inputs, masses, cuts, source intervals and full20ns target remain unchanged.
The report records the requested shift; impossible scales reject. This addresses
measured arithmetic heat loss, not missing physical thermal transport. The
actual-source100digit comparison and complete20ns control are in
`docs/VARIANCE_TRANSFER.md`.

`--coulomb-log` optionally enables the existing displaced-Maxwellian paired
charged-species collision source, with the supplied positive logarithm. The
default remains reaction-only. Each local source half uses the existing
collision-half/reaction-full/collision-half composition; updated species
momenta and energies return to their actual moving supports and alter subsequent
current, Gauss field and pressure. All occupied charged pairs participate once.
No prescribed drag or duplicate heating is applied. Ionization, recombination,
physical masses, exact vacuum, full20ns target and recorded transport steps stay.

This joins previously independent moving/reaction and homogeneous collision
components. The latter's complete independent velocity-integral comparisons are
in `docs/ATOMIC_COLLISIONS.md`; they do not validate the new coupled trajectory.
The declared logarithm, Maxwellian shape, small-drift atomic source, missing
neutral channels/screening and current spatial accuracy remain limitations.
Use the matched full reaction-only trajectory to measure the collision effect;
the low-density slab is not the provisional chamber-fill experiment.

`--neutral-elastic-file PERSONAL_LXCAT_TEXT` additionally enables sourced Ar
ground-neutral/electron elastic exchange, requiring explicit `--coulomb-log`
for this combined source. Exactly one Ar ELASTIC momentum-transfer table is
selected. Both species receive the same paired force/heat with recoil; their
physical masses and actual populations remain. Other neutral collision channels
are still absent. The report records private input SHA, attribution, energy bounds
and maximum observed probability outside the table at collision stages; that
probability is not a missing-force/heat bound. No provider table is bundled.

`--density-jump-minimum-contrast FACTOR` enables existing-face velocity
refinement at represented density jumps after each source stage. The optional
`--density-jump-species ID ...` limits it to those atomic species; omission
selects all represented species, including electrons. Use explicit heavy
indices0..18 to isolate heavy refinement; `--refine-heavy-species` controls
midpoint enrichment only. Contrast must be finite>1 and is a spatial control,
not a density cutoff. Every old knot, fine cell, species and full20ns target
remains. Reports preserve both settings and each block's actual additions.
Defaults remain disabled. The full Ar9-only5ns-source evaluation is saved in
reactive_material_argon_density_jumps_5ns; its remaining Ar8 error motivates
the separate all-heavy experiment, not a general accuracy claim.
Full all-heavy10ns-source contrasts10000/100 now complete20ns in1024steps;
their outputs and measured limits are in all_density_jumps_10ns and
all_density_jumps_c100_10ns. The captured5ns stage remains a separate live run.
