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
