# Full shock and plasma evolutions with explicit compression refinement

Evaluate the concrete cell-bisection option in compressive_refinement.md using
the existing compatible MaterialSpeciesSet equations and minmod numerical
viscosity(.5,2). `--case shock` uses the established finite free argon Sod slab:
half-width50mm, left/right n=1e19/1.25e18m^-3, T=1/.8eV, full duration
0.12*a/sqrt(eV/m). Compare |x|<a/2 with independent gamma5/3 Euler Riemann
density, velocity and pressure, before exterior rarefactions enter that region.

`--case plasma` uses the original physical-mass ion/electron problem, n=1e14m^-3,
half-width10mm, Ti/Te=.2/2eV, zero initial velocity/field and full20ns target.
Cosine initial labels retain free-edge resolution. Only electrons refine;
shock refinement selects its neutral material. The explicit plasma-phase cap
is fixed at .01 by default, independent of subsequent mesh growth. All actual
acoustic limits still apply. No reaction, collision, wall or formation claim.

Before each transport step, bisect all selected cells until the compression
criterion is satisfied or the explicit child-width limit is reached. Keep all
old material nodes. Requested capacity overflow, negative heat or invalid
geometry fails the original full-duration case with its accepted endpoint.
No quiet retry, shortened success, population cutoff, mass change or repair.

Save initial/final complete states, every actual accepted dt and energy/number/
momentum history, and each refinement event's counts/conservation measurements.
This retains deterministic full replay from initial data, source and actual
steps without storing a padded variable-mesh state at every timestep. Failed
cases preserve their accepted prefix and terminal failure. The reported error
norms use the actual requested endpoint only on success. No unit tests.

Usage: `python -m examples.compressive_refinement FRESH_OUTPUT --case shock
--cells 64 --threshold .02 --minimum-child-width 0.00005`. The width is a
declared resolution limit, not a cell-width floor applied to the state.

The first shock attempt saved its endpoint and complete step/budget arrays but
failed JSON export because the success flag was a NumPy boolean. The flag is
now explicitly converted to bool; no dynamics change. Preserve that original
failed attempt and use a fresh full-duration replay for complete event records.

`--thermal-reconstruction limited_linear` selects the explicit bounded pressure
reconstruction documented with the helper. Default constant preserves the
original experiment. Fresh output directories preserve both outcomes. Record
the selector and each event's thermal redistribution/temperature change along
with conserved totals. Evaluate full neutral shock and plasma trajectories
before accepting this candidate; the initial constant-child method failed
the current-accuracy objective despite smaller electric-field differences.

`--no-refinement` disables the refinement loop entirely for a full fixed-label
reference. The report records refinement_enabled=false; threshold, width and
thermal-reconstruction settings are then inactive metadata. Material faces still
move with the same physical velocity equations. This permits the next2048-cell
reference without padded full histories; initial/final states and every actual
dt and conserved budget remain. Use the original phase rule .04*32/N, giving
.000625 at2048. This continues a numerical convergence comparison and does not
turn the finest available trajectory into an exact physical baseline.
