# Moving reaction transfer on precise supports

`PreciseReactionOverlap` keeps the inherited actual atomic-moment evolution and
physical observations while gathering on decimal-coordinate union intervals.
Points and interval membership are resolved before float physical coefficients
are formed. Empty exterior rows remain exactly zero. Physical masses, charges,
donor moments, initially empty products and source normalization retain their
existing checks; no population is discarded to simplify the geometry.

`project_precise_reactions` integrates all original quadrature weights into
common piecewise-constant species amounts and energies, solves the consistent
P1 momentum projection, and recovers each actual U=E-K. Negative thermal states
reject. All source fields remain float64 with per-species scaled N/U; precise
positions and nested fractions carry the distinct support geometry.
The mass solve uses `mass_scaling.solve_mass_banded`: diagonal equilibration
of the same physical matrix and momentum load. This removes the demonstrated
rare-species arithmetic failure without altering density, mass, velocity
resolution or thermal recovery. Negative recovered U still rejects.

Only the explicitly requested earlier ULP cleanup is available. Its tolerance
is coalesce_ulps times the float64 spacing at the remote-domain scale, declared
separately from50-digit geometry. Every species occupancy endpoint is protected.
There is no physical minimum-width coarsening option. Cleanup can still alter
represented local fields; old/new Gauss and field energy differences are observed
without compensating heat or subtracting the residual from physical energy.

`PreciseReactiveMaterial` inherits the original source-half / material-full /
source-half and complete budget contracts. It replaces gather/return only.
It preserves actual variable populations and supports throughout motion rather
than keeping a frozen reaction fixture. General charged accuracy, mesh growth,
source dynamic range and all original formation gates remain separate issues.

By default, source return retains every existing velocity knot still inside
the resulting species support, including old endpoints, and adds the new
support endpoints. These knots are protected from optional ULP cleanup. Newly
occupied species/blocks without old knots use the declared initial cosine
resolution. The actual velocity-cell count can grow; no inherited knot is
discarded to meet that initial count. `preserve_velocity_knots=False` retains
the historical rebuild option for explicit comparisons.

This addresses the measured zero-source transfer error: rebuilding cosine
velocity grids changed final current by55–65% in complete20ns controls despite
near conservation. Retaining knots preserves existing velocity-space detail;
it is not a positivity guarantee for newly mixed or newly born material.
Small physical extensions can impose a shorter actual acoustic timestep,
which must be reported rather than bypassed. Local thermal rejection remains.

`last_velocity_knots` records the exact inherited coordinates used by the latest
projection. It and `last_overlap` are cleared before each source stage, so a
failure while gathering/evolving cannot mislabel an older source as rejected.

`energy_projection='local_total'` remains the default. The explicit research
alternative `'projected_frame'` uses projection_energy.md: a local energy
redistribution paired with the projected momentum, conserving global energy
by Galerkin orthogonality. It retains all density and velocity knots and never
clips thermal energy. Numerical thermal gain, local redistribution L1 and
maximum cell-temperature change are recorded. These errors must converge;
this option is not validated formation physics or physical atomic heating.

Handoff evidence (2026-10-16): the explicit candidate succeeds on the saved
retained-knot10ns rejected source, but Ar9+ reaches2.882045e6eV from~0.199eV.
Global energy error3.03e-15 relative therefore does not justify adoption.
See outputs/reactive_material_argon_retained_knots_10ns/projected_frame_local_errors.json
and docs/THREAD_HANDOFF.md. No full candidate trajectory has been run;
local_total remains the default. Beads5.5.4/.5 track the local velocity and
energy defects before full E2E refinement under5.5.6.

`projection_method='bounded_moments'` is a separate explicit research transfer
with paired momentum/energy corrections and positive cell/node redistribution;
see bounded_reaction_projection.md. It uses the same source, density cuts and
velocity knots. It requires the default energy selector because its paired
energy rule replaces both Galerkin energy choices. Per-block limiter activity,
numerical thermal change and local energy redistribution are recorded. The
default `projection_method='galerkin'` and its rejection behavior are unchanged.

`projection_method='constrained_velocity'` minimizes the volume-norm velocity
error with exact block momentum and convex fine-cell kinetic-energy constraints.
It recovers U=E-K without changing any fine cell's total energy. See
constrained_velocity_projection.md for its strict feasible-start limitation,
barrier tolerances and diagnostic records. This is explicit research code;
local cooling toward an active energy constraint is an approximation error to
measure, not evidence of physically correct heat transfer.

Explicit extend_velocity selects inherited_velocity_faces from velocity_extension.py. Existing old nodes are retained; short newly occupied edge extensions continue the old affine basis. Bounded_moments rejects negative extrapolation shapes. Defaults remain unchanged.
# Explicit newborn mesh option

`newborn_spacing_fraction` is passed to new material models. Retained meshes
bypass selection; density cuts and physical endpoints are unchanged. See
`velocity_selection.md` for the declared target-relative rule and reduced
newborn dimension. The default remains zero. Projection records retain actual
per-block velocity-cell counts, so this option cannot silently claim the
requested dimension was achieved.

`projection_method='variance'` pairs an occupied-volume fitted velocity shape
with direct source heat plus squared velocity-mixing error. Its scalar shape
factor enforces block energy while the shape mean enforces momentum. See
`variance_projection.md` for the derivation and limitations. Local energy
redistribution is recorded, never reset from a global residual. All methods
now report maximum temperature change and minimum temperature per block.
Defaults remain unchanged; this alternative requires trajectory evaluation.

`velocity_refinement_rounds` optionally enriches each inherited/newborn mechanical
mesh before solving its moment transfer, using `velocity_refinement.md`.
All old knots and density cells remain; actual additions and capacity-limited
rounds are recorded. `maximum_refined_velocity_cells` bounds additions by whole
rounds and never shrinks an inherited space. Default zero rounds leaves existing
behavior unchanged. The acoustic bound uses the actual enriched mesh.
Optional `velocity_refinement_species` restricts additions to declared atomic
species indices; all other species still evolve and retain all old knots/cuts.
None selects all species. This separates heavy-species thermal resolution from
electron transport resolution; it never removes a species or changes its mass.

# Explicit binary coordinate scaling

`normalization_shift_bits=b` stores N_relative, momentum and thermal moments
multiplied by2**b and number_scale divided by2**b. This is the same physical
state and same projection equations; the default b=0 is unchanged. Apply the
exact binary shift after dimensionless species normalization and before SI
constants or fine quadrature weights. This preserves positive mixing heat that
otherwise underflows for retained subnormal cell amounts. Requested shifts
are integer0..512; an unrepresentable physical scale rejects instead of clipping
a population or substituting a different shift. Per-block records expose b.
A single scale still has finite dynamic range; this is not arbitrary precision
or a proof of accuracy for the source chemistry. Actual-source100digit moment
integrals and complete trajectories evaluate the change.
