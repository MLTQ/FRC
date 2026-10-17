# Zero-source moving transfer control

Runs the original full20ns planar neutral/ion/electron initial condition at32
velocity cells, physical argon/electron masses, densities1e14/m³, half-width
10mm and temperatures0.1/0.2/2eV. DI/RR are explicitly disabled in every case;
this is an independent zero-source limit, not a replacement reacting model.

Direct transport at half the actual acoustic/plasma bound supplies the control.
Two otherwise identical calculations insert the existing reaction gather/return
every10ns or5ns with a zero reaction source. All distinct density cuts, exact
vacuum and the full original endpoint remain. Changes then arise from the
numerical transfer/changed discrete spaces and step grouping, not reactions.
Energy conservation alone does not establish that this transfer is accurate.

The runner saves initial/final states, actual internal steps and physical
energy, source archive and terminal results. Equal-time electric/current L2
differences use two-point Gauss integration on the exact union of all faces.
An incomplete case is reported as failed and is not compared at unequal time.
This focused end-to-end control diagnoses the actual charged solver change;
it adds no unit-test infrastructure. No coarsening, clipping, thermal repair,
adjusted masses or formation claim.

The initial rebuilt-cosine run is retained in outputs/reaction_remap_control.
Subsequent runs use the current production transfer, which inherits existing
velocity knots. A failed direct control leaves no reference; another remapped
case is never silently substituted as the direct baseline.

`--projection-method bounded_moments` evaluates the explicit bounded paired
moment/energy transfer against the same full direct trajectory. Both numerical
thermal change and energy redistribution are saved alongside J/E differences.
Default Galerkin behavior is unchanged. This comparison measures any diffusion
introduced by the candidate even when the atomic source is exactly zero.

`--projection-method constrained_velocity` selects the volume-norm velocity
projection with physical block momentum and each fine cell total energy
constrained. Source and transport physics are unchanged; barrier/constraint
activity is a numerical approximation, not an adoption criterion.
