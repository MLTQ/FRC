# Combined spatial refinement and complete trajectory replay

The full FRC experiment-sandbox goal remains open. This continuation uses only
end-to-end verification: no unit tests were created or run. See
JOINT_E2E_DESIGN.md for the unchanged equations, physical assumptions and
boundaries. The inspected scientific figure is
outputs/joint_e2e_report/verification.png, with a vector PDF alongside it.

## Evidence locations

All paths below are relative to the repository root. A terminal rejected
workflow is recorded separately from a successful physical calculation.

| Evidence | Directory |
|---|---|
| Both combined grids and full simulation replay | `outputs/radial16_direct_joint`, `_replay`, `_verification` |
| Finest amplitude and tighter-bound controls | `outputs/radial16_finest_controls` |
| Direct physical current/E amplitude difference | `outputs/radial16_current_amplitude_difference` |
| Signed 100 μs trajectories | `outputs/split_mode_e2e_t60_minus` and `_plus` |
| Full trajectory reruns | Each trajectory's `_replay_v2` and `_verification` |
| Independent trajectory field/conservation audit | `outputs/split_mode_e2e_t60_audit` |
| Signal-energy normalization | `outputs/split_mode_e2e_signal_energy` |
| Scientific figure and plotted inputs | `outputs/joint_e2e_report` |
| Interrupted work, source inventory and resource estimates | `outputs/radial16_joint_diagnostics` |

## Combined spatial response

Every row below uses the same independently resolved radial16,N128 reference,
growth12307.1561445/s, actual signed nonlinear RHS and physical quadrature20.
F denotes fluid radial elements/axial points; B denotes field cells/axial points.
The complete new combined-grid cases and their fresh full reruns agree in
every saved coefficient, physical fingerprint and non-runtime metric.

| F | B | Radial velocity rate | Radial B rate | Toroidal B rate | Current rate |
|---|---|---:|---:|---:|---:|
|32/128|128/128|0.1088365|0.05193797|0.3830690|0.08756878|
|32/128|256/128|0.03090839|0.05162497|0.2793811|0.03375043|
|32/256|256/256|0.03088541|0.000982554|0.04947849|0.02591733|
|64/256|256/256|0.001831426|0.000979624|0.005475095|0.000878993|

Values are relative physical L2 errors, not percentages or coefficient norms.
Fluid refinement on the already refined field grid reduces current-rate error
29.5-fold, radial-velocity error16.86-fold and toroidal-field error9.04-fold.
The finest full electric-field response error is4.568074e-6 and density-rate
error5.160122e-5. Worst strict induction residual is6.461405e-11 under the
unchanged1e-10 bound. No numerical diffusion, filtering or repair was added.

On the refined field grid, initial weak-current error is4.190757e-6, strong-
Ampere discrepancy5.475849e-5 and internal radial normal-current jump1.033910e-5.
Accurate initial fields alone therefore do not certify the much more sensitive
rate response. Strong volume curl still excludes distributional face terms.

The new quadrature16/20 norm differences are at most5.00e-15 absolute. The
largest relative difference10.44% belongs to an initial axis-density error
4.3–4.8e-14; its absolute difference is5.00e-15. This precision-scale outlier
is retained rather than described as universal relative quadrature agreement.
For F32/256, the full Faraday rate-relative diagnostic reaches3.506e-7 while
its original backward error is below3.55e-17. Both diagnostics remain visible.
The spatial quadrature control and its complete replay passed.

## Complete time evolution

A separate radial1,N64 contact benchmark runs both signed perturbations for
100us in60 guarded Gauss steps on F4/32,B8/32. The independent mode supplies
only initialization and audit targets; actual nonlinear fluid, species,
current and magnetic variables evolve. New opt-in split_mode_evolution uses
the split-inverse arithmetic with every original stage, endpoint, positivity,
induction bound and iteration budget unchanged.

Both trajectories and complete fresh reruns finished. Every accepted stage,
state and non-runtime metric matches bitwise:62 NPZ records per sign. Worst
stage induction residual is9.989630703e-12, maximum nonlinear residual/limit
0.9604443, and maximum endpoint relative residual1.63839e-16.

Independent quadrature16/20 reconstructs all fields, pressures and conserved
quantities. Final quadrature20 relative physical errors are:

| Quantity | Error |
|---|---:|
| Density |0.00101070362|
| Radial velocity |0.00203901899|
| Radial magnetic field |0.00193651388|
| Toroidal magnetic field |0.00065585717|
| Current |0.00274790798|
| Electric field |0.00867838813|

The independent growth factor is2.402175837; component projections range
2.401218278–2.402203516. Reconstructed final fields agree exactly with the
older same-grid60-step trajectory in the saved comparison. Both states' E
fields use the new arithmetic, so this is not a claim about identical internal
arithmetic between the two solver sources.

Initial physical perturbation energy is1.652636476e-11J. Actual total-energy
drift divided by that scale is-1.315005e-7 and-3.635322e-7 for minus/plus at
quadrature20; quadrature16 gives-1.314897e-7 and-3.635445e-7. The scale is not
a conserved free energy. Absolute energy changes are-2.173225e-18J and
-6.007865e-18J. Particle relative drift is at most1.75e-17; specific-entropy
drift at most3.52e-17 per particle. Sampled density stays above7.4997e19/m³,
and electron/ion temperatures above3.9648/0.3722eV. The physical audit,
old/new comparison and energy normalization all replay completely.

These coarse integrations do not establish full refined time/mesh convergence.
The previously observed unwanted growing branch and failure to demonstrate
fourth-order time convergence remain open requirements.

## Equivalent computation with lower setup and audit cost

All earlier Python sources remain unchanged. Each improvement is opt-in.
Complete signed comparisons isolate each change; the adopted variants also
have complete fresh replays:

- Radial slicing retains every quadrature point and physical fingerprint.
  All raw arrays and solver diagnostics match the older full-array baseline;
  changed norm reduction order gives at most2.22e-16 absolute difference.
- Reusing the identical axial curl basis preserves every metric and fingerprint.
  Observed baseline runtime decreases245.685s to98.613s.
- DirectMetricWallPlasma skips only the magnetic-Hessian LU that the original
  subclass immediately discards. Coordinates, used metric and full solve code
  remain identical. The complete large F32/256 case agrees in every raw array
  and diagnostic; observed runtime decreases707.845s to480.598s.
- A within-run cache reuses only original immutable field geometry. Every
  actual fluid space, state, density, species and induction operator remains
  fresh. Full original/new controls with two fluid grids sharing one field
  grid agree bitwise, and the entire cached batch replays.

Runtime comparisons are observations, not controlled performance benchmarks.
None changes a factorization algorithm used by the actual solve, a physical
mass, accepted state, boundary, residual bound or diffusion model.

The first full-array large audit was interrupted during reconstruction when
the host reported no writable disk space. A concurrent first trajectory replay
exited1 after two saved steps without a logged cause. Later setup transitions
also left explicit partial batches. Every interrupted archive remains intact.
Completed-case selections record their parent batch's incomplete status and
verify only the unchanged selected case, never unfinished siblings. Fresh full
trajectory/spatial runs subsequently completed. See interrupted_workers.json.

## Finest amplitude and residual controls

The F64/256, B256/256 smaller-amplitude case has completed at epsilon 1e-5
under the original 1e-10 bound. Its maximum induction residual is 6.935688e-11.
Quadrature 20 errors compare as follows:

| Rate | epsilon 1e-4 | epsilon 1e-5 |
|---|---:|---:|
| Radial velocity | 0.00183142646 | 0.00183390466 |
| Radial magnetic field | 0.00097962351 | 0.00097962352 |
| Toroidal magnetic field | 0.00547509477 | 0.00550282935 |
| Current | 0.00087899301 | 0.00088543967 |
| Electric field | 0.000004568074 | 0.000004568896 |

Axial-velocity error changes from 7.53577e-7 to 9.22389e-7, so relative
agreement is not uniform across components. Error-norm changes alone do not
measure the physical norm of the difference between the two rate fields.
Centered amplitude differences in Euclidean coefficient norms are
6.29740e-5 for the fluid rate, 1.76520e-4 for magnetic rate, 0.00947096 for
current rate, and 6.58187e-8 for E. In particular, the 0.947% current coefficient
change remains visible despite the close physical error norms. Coefficient
amplitude convergence is not established by these error-norm comparisons.

The additional `current_amplitude_audit` reconstructs both accepted responses
independently in physical space. Every current/E/reference/weight fingerprint
matches the original audits; no plasma solve or field-mass factorization is
needed. At quadrature 20, actual field-difference L2 norms divided by the
independent reference norm are:

| Centered response | Relative physical difference |
|---|---:|
| Full current-rate vector | 0.000107450107 |
| Radial current rate, relative to its own component | 0.0301482272 |
| Azimuthal current rate, relative to its own component | 4.528043e-10 |
| Axial current rate, relative to its own component | 0.00126286432 |
| Full electric field | 7.148893e-8 |

Thus a 0.010745% total current difference hides a 3.0148% change in the much
smaller radial component; the axial component changes 0.1263%. These are actual
physical differences, distinct from the 0.947% coefficient change and from
changes in error norms. The two quadrature orders' relative L2 differences
agree to within 8.68e-19 absolute. Sampled maxima are retained separately and
vary with quadrature points; no continuous maximum bound is claimed. Absolute
values refer to centered unit-amplitude responses. The complete reconstruction
replay matches every saved array, fingerprint and metric bitwise.

The epsilon 1e-4 case with a tighter 1e-11 bound rejects the minus input before
returning any rate. Across five inner bounds from 1e-11 through 1e-15, the final
dual Faraday residual stays between 5.8882e-11 and 5.8939e-11. Operator residuals
remain 3.2857e-11 to 3.3944e-11 and physical Faraday residuals about 1.588e-11;
Ohm and Ampere backward errors remain below 5.39e-16. Every original refinement
budget is retained. This is a measured tighter-bound rejection, not a completed
physical comparison or proof of its floating-point cause. Its exact saved input
has SHA ce34210139b000f4c2746dfde664ba1dd188d5e2485572a364e184bf12c4d999.
Every saved state and initialization array is bitwise identical to the accepted
same-grid input at the original bound; no failed rate array was returned.

The two residuals share the same right-hand-side normalization. With
`u = operator(solution) - rhs` and
`v = dual_lift.T @ (faraday @ electric - rate)`, exact arithmetic gives `u = -v`.
The reverse triangle inequality bounds their normalized mismatch from below by
the absolute difference of their reported norms: 2.499e-11 to 2.607e-11 at the
recorded final iterates. Each lower bound exceeds twice the requested 1e-11
bound. This identifies incompatible evaluated residual paths at those iterates;
it does not establish impossibility at another iterate or isolate the cause.
The saved derivation is `tighter_residual_bounds.json` in the joint diagnostics.
Next capture the actual vectors and individual split-inverse actions on this
exact input. The old balanced-arithmetic observation probe is not a substitute
for an observation-only copy of the current split-inverse arithmetic.

The original batch and fresh complete rerun are terminal with one complete and
one failed case each. Every saved array, failure message and non-runtime metric
replays bitwise. The paired audit and its replay explicitly mark the tighter
case `not_comparable_failed_solve`. The 256 MiB disk-headroom guard did not
trigger. The direct physical current/E difference audit above measures only
those fields; fluid/magnetic physical difference norms remain separate work.

## Remaining physical scope

The filled periodic argon fixture is not a selected laboratory experiment.
The user's approximate1m quartz tube,10–20cm diameter, selectable gas,
0.1–100mTorr neutral fill and adjustable bias remain provisional. Full refined
trajectory convergence, local budgets, heated moving exact vacuum,
transport/ionization, real ends/coils, converged formation and experiment UI
remain unfinished. Axisymmetric closed flux establishes neither3-D stability
nor particle confinement. All 12 acceptance beads and the overarching goal stay
open. No remote or publication is configured or used; the user's poster is
preserved unchanged.

## Final checkpoint

All workers are terminal and collected. Twenty complete verification groups
cover simulations, exact rejection replay, physical postprocessing and an
unchanged completed-case selection; this is not a count of successful physical
experiments. Their 137 raw-file comparisons and non-runtime metrics match.
All 71 new source archives, 217 preceding archives and 22 earlier source closures
are verified. Fourteen new modules each have companion documentation; no older
Python file changed. Final source union: 128 files, SHA
8a97cac82e67e4aa62eeea9f5c80edfd1c982529a158edd6fac0512d9f56e457.
See `outputs/radial16_joint_diagnostics/checkpoint.json` for the exact inventory.
The user's poster and all interrupted evidence remain unchanged. No unit tests
were added or run, and no remote or publication was used.
