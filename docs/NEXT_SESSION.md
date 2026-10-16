# Resume the physical FRC simulation

## Paired pulse time sensitivity and published heat coefficient (2026-10-16)

The previous user turn was a UI recheck with no new scientific progress.
This continuation collected the completed paired half-timestep pulse and used
its field differences to launch a necessary third full timestep level. The
original sandbox objective and all remaining acceptance gates stay open.

Paired run33952 completed8000steps/3.5us (5172.282s,exit0), with all8000 raw
states and matching final digest. Source remains31721f1fef01b747862f5263f3ddb92e5dfbca0d3fc1995e5e03ebf229d17afc.
Existing actual-basis comparison13672 exited0; PNG inspected, PDF saved in
outputs/material_tube_maxwellian_timestep_comparison. Onlydt differs from the
4000-step run. Final vector differences: B10.1933%, J20.9077%, velocity2.4197%,
maximum position mismatch4.32214mm. Sampled J difference peaks62.3226% at2.58125us.
Collision heat changes109.0123 to107.8718mJ; numerical ion heat32.0043 to32.3384mJ,
far above fine physical ion collision heat0.43423mJ. Max work residual/E0 improves
6.12049e-4 to8.95750e-5, but detailed fields are NOT time-converged. Read
FINITE_DRIFT_COLLISIONS.md. No unit tests or new framework.

New paired16000 configuration is committed as60f76e4e4 on the isolated branch
codex/finite-drift-collisions at /private/tmp/frc-finite-drift. The full3.5us
pulse is running at0.21875ns with unchanged physical and numerical-pressure
inputs, source31721f1f...d17afc. All results stay in the main outputs directory.
Main pulse source remains05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75.
Do not merge the paired moving code while the original main jobs remain live.

The completed kinetic operators independently reproduce the published Z=1
stationary leading nonlocal coefficient b=43.5: degrees20/24/28 give43.59610 /
43.56493 /43.54514. Finest difference0.1038%, no fitted matrix. One-off analysis
and input hashes are in outputs/kinetic_heat_published_comparison. Read
KINETIC_HEAT_TRANSPORT.md for normalization and scope. The paper's short-wave
SPRING curve is a complex modal conductivity modulus, not our instantaneous
real transient heat/temperature ratio; no short-wave published accuracy claim.
This remains a uniform periodic limit, not heat transport through a plasma edge.

FOUR pulse handles were actually polled LIVE. Most recent saved observations:
-16339 paired4x16t16000:400/16000,257.895s;0.21875ns,
 outputs/material_tube_independent_20mt_maxwellian_r4_z16_t16000.
-49353 paired6x24t8000:2000/8000,3102.810s;0.4375ns,
 outputs/material_tube_independent_20mt_maxwellian_r6_z24_t8000.
Both use worktree source31721f1f...d17afc.
-58668 original4x16t16000:4200/16000,2625.143s;0.21875ns,
 outputs/material_tube_independent_20mt_viscous_r4_z16_t16000.
-85554 original6x24t8000:5100/8000,7948.597s;0.4375ns,
 outputs/material_tube_independent_20mt_viscous_r6_z24_t8000.
Both use main source05555dc9...b3c75. Saved counts are observations, not liveness
proof; re-poll exact handles. Never restart from an observation timeout.

NEXT collect complete original/paired16000 and spatial8000 results and use the
existing --configured-fields comparison from the isolated worktree. Keep source
closures intact. In parallel, continue the missing nonuniform heat evolution
with finite current and explicit edge/end energy/particle assumptions; the
periodic zero-current infinite-ion-mass limit cannot be installed into formation.
No converged formation, physical transport, stability or confinement is claimed.
Keepfrc-transport.3/.4,frc-vacuum.3.3 and the full goal open. Existing10ns viewers
and the user poster are unchanged; all research stays local with no remote.

## Linear kinetic heat transport completed (2026-10-16)

Previous goal turn was progress:66ca30ac2 saved complete paired collision and
original half-timestep pulses, exposed large current sensitivity, and launched
an additional full timestep refinement. This continuation kept all four pulse
runs live and their scientific source closures unchanged while evolving the
missing heat-transport physics in an independent uniform-background limit.

New plasma_bench.kinetic_heat retains full linearized electron-electron Landau
collisions, velocity-dependent electron-ion Lorentz scattering and dynamic
velocity moments. The periodic small-temperature-perturbation fixture assumes
stationary infinite-mass ions and a quasineutral zero-current electric response.
It is NOT installed as the moving pulse heat closure. No BGK substitution,
flux limiter, speed/collision floor, artificial damping, state repair or unit
suite. Read KINETIC_HEAT_TRANSPORT.md and both new companion documents.

Complete heat calculations89659 and92079 exited0. The refined output is
outputs/kinetic_heat_relaxation_refined, source
acb28fcfc9834a0be3d29ca5e96eed07747f6a1b255394ee12580a99275e16cc.
Both scientific PNGs inspected; PDFs, complete operators and state histories
saved. Degrees8..28 improve the independent analytic Lorentz conductivity
comparison to2.4604e-6 relative difference. Full Z=1 local conductivity reaches
3.20203 in Braginskii units,1.33% above the approximate3.16 reference, without
fitting. Collisionless histories converge against independent Gaussian-
characteristic Volterra evolution: degree28 max temperature error3.4723e-8,
near the independent reference-integration floor; free-energy drift2.89e-15.
The strongest collisional case K=1 still changes by1.696e-3 of initial
perturbation amplitude between degrees24/28 (about1.65% of remaining final
amplitude), so do not claim all kinetic regimes converged. Nonlocal heat
relaxation is much slower than local diffusion in these complete trajectories.

Code saved on branch as e275ffe57; ONLY its four new standalone files were
cherry-picked to main as469927bc8. Main heat source exactly matches the refined
run; main pulse source remains05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75.
The paired moving collision source31721f1f...d17afc still resides on isolated
codex/finite-drift-collisions at /private/tmp/frc-finite-drift; do not merge that
coupling while original main runs remain live. Existing10ns viewers unchanged.

FOUR exact pulse handles were re-polled LIVE; these are observed saved counts:
-58668 original4x16t16000:1800/16000,1138.54s;0.21875ns/full3.5us requested.
-85554 original6x24t8000:4100/8000,6424.51s;0.4375ns spatial control.
Both main source05555dc9...b3c75, in
outputs/material_tube_independent_20mt_viscous_r4_z16_t16000 and
outputs/material_tube_independent_20mt_viscous_r6_z24_t8000 respectively.
-33952 paired4x16t8000:6700/8000,4342.47s;time control.
-49353 paired6x24t8000:1100/8000,1723.27s;spatial control.
Both worktree source31721f1f...d17afc, in
outputs/material_tube_independent_20mt_maxwellian_r4_z16_t8000 and
outputs/material_tube_independent_20mt_maxwellian_r6_z24_t8000 respectively.

NEXT collect those exact handles; never restart from an observation timeout.
From /private/tmp/frc-finite-drift with MPLBACKEND=Agg and writable caches,
use existing --configured-fields comparisons on complete successful pulse pairs,
including paired4000/8000, original8000/16000 and matched8000 spatial grids.
The independent kinetic baseline now provides actual local/nonlocal histories
for the next heat model, but nonuniform/finite-drift/moving heat evolution and
explicit edge/end energy/particle budgets remain. Do not transplant its periodic
zero-current/infinite-ion-mass assumptions into formation. frc-transport.4 is
now in_progress; full goal,frc-transport.3 andfrc-vacuum.3.3 remain open. Preserve
poster and live sources, finish numerical comparisons, then integrate defensible
transport with the full moving physical equations. No remote/publication.

## Completed paired pulse and unresolved timestep sensitivity (2026-10-16)

The preceding UI-only turn made no new scientific progress: it rechecked an
already delivered10ns replay. This continuation collected complete physical
results and changed the next numerical action from evidence. The full original
objective and acceptance beads remain open; no unit suite or new framework.

Paired-collision4000-step run52937 completed3.5us (2621.768s,exit0). All raw
states and final digest are present. Existing physical comparison38724 and
actual transport review90288 completed, with all three scientific PNGs inspected
and PDFs saved. Read FINITE_DRIFT_COLLISIONS.md for the complete measured table.
At final time paired collision heat is109.0123mJ versus238.0457mJ prescribed eta;
particle RMS radius is22.9773 versus25.3185mm. Both expand from18.7082mm.
Numerical ion heat is32.0043mJ versus only0.43347mJ actual paired ion collision
heat; this is not a physically validated ion-temperature prediction. The paired
max observed work residual/E0 is6.1205e-4; its time comparison remains pending.
At peak drive47.420% of particles and81.247% of actual drag power are in the
conditional electron parallel density-Kn review region. Spatial/nonlocal heat
transport remains required; fixed-log Maxwellian drag is not full transport.

Code for actual/hypothetical collision labeling, physical comparisons and the
new paired6x24 configuration is committed as f5a7bb66d on local branch
codex/finite-drift-collisions at /private/tmp/frc-finite-drift. Scientific runner
source remains31721f1fef01b747862f5263f3ddb92e5dfbca0d3fc1995e5e03ebf229d17afc.
Main-checkout original solver remains unchanged. All outputs stay in the main
/Users/max/Code/FRC/outputs directory; no remote exists or publication is needed.

Original half-timestep run88009 completed all8000steps/3.5us (5189.560s,exit0).
All raw states/final digest present. Existing configured field comparison18763
and independent peak-linkage calculation88092 completed(exit0); PNG inspected,
PDF saved in outputs/material_tube_array_timestep_comparison. ONLY dt changed
from0.875ns to0.4375ns. Across81 common times/6400 material samples, final vector
L2 differences are B8.240%,J20.583%,velocity3.030%,position mismatch2.463mm.
Current difference reaches71.988% at1.575us despite similar heat/RMS radius and
small energy errors. Detailed fields are not time-converged. Read
MATERIAL_COMPRESSION.md. A full third level at0.21875ns/16000steps was therefore
launched, preserving all physical/numerical terms and full3.5us duration.

FOUR exact scientific handles were polled LIVE; counts below are observations:
-58668 outputs/material_tube_independent_20mt_viscous_r4_z16_t16000,
 first accepted step;0.21875ns,16000steps/full3.5us requested.
-85554 outputs/material_tube_independent_20mt_viscous_r6_z24_t8000,
 3400/8000steps,5342.11s;original spatial control.
Both use main source05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75.
-33952 outputs/material_tube_independent_20mt_maxwellian_r4_z16_t8000,
 4900/8000steps,3202.66s;paired time control0.4375ns.
-49353 outputs/material_tube_independent_20mt_maxwellian_r6_z24_t8000,
 300/8000steps,470.28s;paired spatial control,6x24/6vacuumcells12modes.
Both use unchanged worktree source31721f1f...d17afc above.

NEXT collect exact handles without restarting on observation timeouts. Use the
existing worktree tool examples.material_coil_plot --configured-fields OUTPUT
FIRST SECOND, with MPLBACKEND=Agg and writable caches, on terminal successful
pairs: original8000/16000 for time; original4x16/6x24t8000 for space; paired
4000/8000 for time; paired4x16/6x24t8000 for space. Do not claim convergence from
two grids. Keep independent linkage, heat and conservation observations separate
from field accuracy. Inspect actual failures if any. Neither checkout's live
scientific source may change. Merge branch into main only after BOTH main jobs
85554/58668 finish, not merely after the now-complete88009. Paired jobs stay in
the worktree. Keep frc-transport.3,frc-vacuum.3.3 and frc-goal open. Pursue
frc-transport.4's physically derived spatial heat evolution after complete paired
comparisons. Preserve the user poster and existing10ns viewers.

## Bulk transport-scale evidence (2026-10-16)

Previous goal turn was progress: main69a72827e saved the complete bias comparison,
raw control trajectory and measured numerical field sensitivity; local branch
350f14014 added the configured actual-field comparison. This continuation kept
both scientific source closures unchanged while completing an actual-plasma
transport review on the original complete long-array trajectory.

Existing material_collision_diagnostic now differentiates the actual n0/det(F)
density map and measures conditional thermal mean-free-path content relative to
density scales, including the direction along B. Full10579 exits0; independent
actual-map derivative comparison83409 exits0, with second-order central-difference
agreement (2.8e-6 ->7e-7 when interval halves). Code50bdadad0 is saved on isolated
codex/finite-drift-collisions, /private/tmp/frc-finite-drift. Data, frozen source,
PNG/PDF and derivative recipe are in outputs/material_array_collision_regime_content;
source68546fde89b65e5a2e705aba9ba766b281888318c794c9b841420384a76998da.
Scientific region-content PNG inspected. Read FINITE_DRIFT_COLLISIONS.md.

At original peak drive1.26875us,47.120% of particles and78.861% of hypothetical
log10 drag power occupy magnetized regions with parallel density Kn>=.1.
Initial weak-field electron density Kn>=.1 contains99.179% of particles; final
parallel review contains92.720%. This is bulk transport-scale evidence, not just
an edge maximum. Conditional NRL thermal rates and density-only gradients do not
prove kinetic validity or predict heat escape. They omit temperature/field scales
and intercell jumps. Newfrc-transport.4 tracks actual directional/nonlocal thermal
evolution with a physical vacuum-edge/end energy treatment. No scalar diffusion,
flux limiter, coefficient floor, atmosphere or altered live trajectory was added.

FOUR exact handles were polled LIVE; saved counts remain observations:
-52937 outputs/material_tube_independent_20mt_maxwellian_r4_z16_t4000,
 2400/4000steps,1626.32s; at2.1us work residual/E0=-6.1201e-4,
 meanTe/Ti12.8714/3.1059eV,collision heat94.2076mJ,numerical ion heat29.1075mJ.
 This is an unfinished pulse; its larger energy residual needs complete time
 comparison, not a correction or truncated-success claim.
-33952 outputs/material_tube_independent_20mt_maxwellian_r4_z16_t8000,
 2200/8000steps,1492.57s;0.4375ns/full3.5us requested.
Both use unchanged worktree source31721f1fef01b747862f5263f3ddb92e5dfbca0d3fc1995e5e03ebf229d17afc.
-88009 outputs/material_tube_independent_20mt_viscous_r4_z16_t8000,
 5600/8000steps,3700.30s; original time refinement.
-85554 outputs/material_tube_independent_20mt_viscous_r6_z24_t8000,
 2300/8000steps,3674.22s; original spatial refinement.
Both retain source05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75.

NEXT collect complete results from those four handles. Do not restart on timeout
or edit their scientific source. Use existing configured field comparison for
original time then spatial controls, and paired4000/8000. Inspect actual failures
if any. Repeat physical region analysis on completed paired states (label actual
versus hypothetical collision power correctly), then pursue frc-transport.4.
Merge local branch into main only after both original-source jobs finish. Keep
full original goal and frc-transport.3/frc-vacuum.3.3 open. No unit tests or
recursive verification framework. Preserve user poster and10ns viewers; no remote.

## Completed bias comparison and measured field sensitivity (2026-10-16)

Previous goal turn was progress: paired collision coupling was committed in the
isolated branch (279bf0a8f), with two full pulses launched and main checkpoint
317818b0e. This turn collected bias-only handle46579 terminal success: all4000
steps /3.5us,2539.966s,source05555dc9...b3c75. All4000 raw states and matching
final digest are present. Comparison51962 exited0; PNG inspected and PDF saved:
outputs/material_tube_array_vs_bias_viscous_comparison. The array final particle
RMS radius is25.3185 versus27.4959mm bias-only (7.919% smaller), but both expand
from18.7082mm. The array final surface is larger and electrons much hotter.
Read MATERIAL_COMPRESSION.md. This completes the matching bias comparison;
frc-vacuum.3.2 is closed, while physical accuracy and the full goal remain open.

Existing material_coil_plot now supports configured pulses at81 common actual
times with higher common quadrature, component B/J/velocity differences and
material-position mismatch. Code350f14014 stays on codex/finite-drift-collisions
in /private/tmp/frc-finite-drift. Complete old/new two-loop comparison92245 exits0
in outputs/material_tube_two_loop_numerical_sensitivity_agg (PNG inspected, PDF
saved). Changed BOTH dt1.75->.875ns and numerical viscosityoff->C1=.5/C2=2:
final differences B14.958%,J25.224%,velocity10.682%,position4.458mm. These are
numerical-treatment sensitivity, not isolated time convergence. Full budget
errors remain small; fields still require accuracy evidence. First plotting
attempt96694 exited134 without results and is preserved with failure.json.
No simulation was restarted or altered because of that postprocessing failure.

FOUR actual scientific handles were re-polled LIVE (saved observations only):
-52937 outputs/material_tube_independent_20mt_maxwellian_r4_z16_t4000,
 1150/4000steps,794.93s; .875ns, full3.5us requested.
-33952 outputs/material_tube_independent_20mt_maxwellian_r4_z16_t8000,
 1000/8000steps,693.59s; .4375ns, matching full pulse.
Both use isolated worktree source31721f1fef01b747862f5263f3ddb92e5dfbca0d3fc1995e5e03ebf229d17afc,
fixed declared log10, paired drag/heat, original4x16 mesh andC1=.5/C2=2.
-88009 outputs/material_tube_independent_20mt_viscous_r4_z16_t8000,
 4400/8000steps,2927.61s; original time refinement.
-85554 outputs/material_tube_independent_20mt_viscous_r6_z24_t8000,
 1800/8000steps,2891.36s; original spatial refinement.
Both retain main source05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75,
verified unchanged. Do not change either checkout's live source or restart from
a timeout. The postprocessing change does not belong to either runner closure.

NEXT collect those exact four handles. Use the existing comparison extension
from the worktree: python -m examples.material_coil_plot --configured-fields
OUTPUT FIRST SECOND, with MPLBACKEND=Agg and writable cache directories. Require
terminal successful inputs; use original4000/8000 for time, then original8000
4x16/6x24 for spatial sensitivity. Compare paired4000/8000 after both complete.
Failures need their actual saved-state diagnosis, not altered physics or clipping.
Merge the local branch into main only after BOTH remaining original-source jobs
finish. All outputs stay in /Users/max/Code/FRC/outputs. No unit suite. Keep
frc-vacuum.3.3,frc-transport.3 and full original goal open. Preserve user poster
and existing10ns /20mt/ and /array20mt/ views; no remote/publication.

## Paired collisions in the full moving solver (2026-10-16)

Previous UI-only turn made no new scientific progress: it verified the already
completed 10 ns animation. This turn resumed the next available physical change.
Local commit 279bf0a8f on codex/finite-drift-collisions adds one paired collision
force/heat calculation to both driven and undriven moving electromagnetic rows
and species entropy. Read FINITE_DRIFT_COLLISIONS.md and the new companion docs.
No old physical terms or exact vacuum were removed; no heat/field/state repair.

The implementation lives in sparse worktree /private/tmp/frc-finite-drift while
three original studies keep their main-checkout source unchanged. The local Git
branch preserves the code; all new output is in /Users/max/Code/FRC/outputs.
No remote exists or is required. Do not edit scientific source in either checkout
while its respective live runs remain active. Main can receive the branch once
its THREE original jobs have terminal evidence; worktree jobs stay independent.

FIVE actual scientific handles were polled LIVE (counts are saved observations):
- 52937: outputs/material_tube_independent_20mt_maxwellian_r4_z16_t4000,
  full 3.5 us requested, 0.875 ns, recorded 150/4000 steps, 102.49 s elapsed.
- 33952: outputs/material_tube_independent_20mt_maxwellian_r4_z16_t8000,
  matching half-timestep pulse, 0.4375 ns, first accepted step recorded.
These two use worktree source
31721f1fef01b747862f5263f3ddb92e5dfbca0d3fc1995e5e03ebf229d17afc,
paired displaced-Maxwellian collisions with declared lnLambda=10, eta=0 and no
separate exchange, original 4x16 mesh and C1=.5/C2=2 numerical bulk viscosity.
Species collision heat and numerical ion heat are independently observed.
At early 0.13125 us, coarse electron/ion exchange is -21.732/+21.732 uJ and
external-work residual/E0=6.67e-15. This precedes pulse onset, not full accuracy.

- 88009: outputs/material_tube_independent_20mt_viscous_r4_z16_t8000,
  original time control, recorded 3400/8000 steps, 2240.10 s elapsed.
- 85554: outputs/material_tube_independent_20mt_viscous_r6_z24_t8000,
  original spatial control, recorded 1400/8000 steps, 2222.76 s elapsed.
- 46579: outputs/material_tube_bias_only_viscous_r4_z16_t4000,
  original matching bias-only control, recorded 3100/4000 steps, 1944.52 s elapsed.
Their main source remains
05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75,
verified unchanged after the isolated implementation.

NEXT collect complete results from those exact handles. Use the existing
material_experiment_compare on completed bias/array controls. Measure complete
field/current/geometry/heat time and mesh differences, not just final maxima.
For paired collisions, compare full pulse and half-dt trajectories and then mesh
sensitivity, retaining declared-log/Maxwellian limitations. Do not restart a live
job because an observation times out. Source baselines remain independent and
complete; these moving pulses remain unfinished. Keep frc-transport.3,
frc-vacuum.3.2/.3.3 and the full original goal open. No unit suite or new test
framework. Preserve the user poster and existing /20mt/ and /array20mt/ views.

## Finite-drift collision model and measured sensitivity (2026-10-16)

Previous goal turn was progress, committed as b4bb4a1da: complete matched coil
comparison, array display and three physical controls launched. This continuation
kept their scientific source unchanged and implemented a paired physical collision
component. Read FINITE_DRIFT_COLLISIONS.md and the new module companion.

plasma_bench/maxwellian_pair_collisions.py evaluates Landau momentum and species
heat moments for displaced isotropic Maxwellians, with explicit positive Coulomb
logarithm. Force is paired and the species heat sum restores lost drift energy;
no drag/heat clipping, floor or atmosphere. This is not yet installed in the moving
coil solver; anisotropy, screening/cutoff and kinetic validity remain unresolved.

Complete source-relaxation run 53955 exited 0. Two cases (equal masses and physical
electron/Ar) at 100/200/400 steps agree with independent velocity-space-integral
trajectories, with roughly fourth-order time convergence. Fine relative-drift
errors 1.390e-7 / 1.285e-7, energy drifts 6.418e-8 / 5.964e-8. Equal masses reach
the independently predicted 5 eV equilibrium; Ar temperatures equilibrate much
more slowly than current. Evidence: outputs/maxwellian_pair_relaxation, source
79b7bf4c446736318c5f970ca50bcce508a9bf269ed27bdd02477f3d4e889080.

Actual-pulse postprocessing 61137 also exited 0:
outputs/material_array_collision_sensitivity, source
ae11182af08eb9de9989fd250a5113268ee9ca50d4c75396e7d707fdc53c6810.
On unchanged original states, declared lnLambda=10 gives finite-drift drag path
integral 63.449 mJ versus 77.678 mJ small-drift Maxwellian and 236.381 mJ original
prescribed eta on the same sparse samples. This is not a changed-transport pulse
or prediction of its eventual heating. Original every-step Joule heat remains
238.046 mJ. Candidate max local relaxation rate times original dt is 0.07222,
which is not a full-system stability certificate. Both new scientific PNGs were
inspected; PDFs generated. No unit tests or new verification framework.

The same THREE actual physical handles were re-polled LIVE:
- 88009: outputs/material_tube_independent_20mt_viscous_r4_z16_t8000,
  recorded 2100/8000 steps, 1315.82 s, nonterminal.
- 85554: outputs/material_tube_independent_20mt_viscous_r6_z24_t8000,
  recorded 800/8000 steps, 1208.26 s, nonterminal.
- 46579: outputs/material_tube_bias_only_viscous_r4_z16_t4000,
  recorded 1750/4000 steps, 1040.76 s, nonterminal.
All keep source 05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75;
current recursive source identity was checked unchanged. Counts are observations,
not completion. Do not restart from timeouts or modify their source while live.

NEXT collect complete physical controls/refinements and compare as below. New
frc-transport.3 tracks replacing prescribed resistance and small-drift exchange
with one paired collision closure in the moving solver. Add/reuse a collision
hook in MaterialVectorCurrent and MaterialCoilCurrent only after live source
restrictions end; avoid double counting or a copied divergent solver. Then evolve
complete physical trajectories with conservation and numerical sensitivity.
The full original goal remains open. Preserve user poster, old displays and all
raw states; no remote publication. Existing /20mt/ and /array20mt/ views remain.

## Completed coil comparison and live physical refinements (2026-10-16)

The previous goal turn made progress: the requested 10 ns animation was completed
and committed as cabe7f572. This continuation collected matched control 84694
(exit 0, 4000 steps / 3.5 us, 2362.285 s). Both physical runs use frozen source
05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75.

The completed comparison is outputs/material_tube_array_vs_two_loop_viscous_comparison.
Read MATERIAL_COMPRESSION.md for the measured table. The array gives 2.776 times
Joule heat but only 1.40% smaller final particle RMS radius than the two-loop case;
both expand overall. Numerical ion heat is separate (15.634 versus 2.691 mJ).
Independent peak-drive coil-flux/current-linkage differences are 0.62–1.26% for
the array groups. At peak drive, drift above electron thermal speed accounts for
7.985% of the array Joule power; it cannot be dismissed from its final edge maximum.
No physical transport, formation, confinement or spatial convergence claim follows.

Completed array display: outputs/material_experiment_independent_20mt_10ns_display,
http://127.0.0.1:8790/array20mt/ through server 90454. Builder 53068 and comparison
3803 both exited 0. There are 351 actual frames at nearest 10 ns targets; 0.875 ns
integration gives 9.625/10.5 ns gaps. Separate bias/array labels and numerical
compression assumptions are visible. Playback, signed current, scales, final
state and disclosure were exercised in the browser with no console warnings or
errors. Scientific PNGs and full-page browser_view.jpg were inspected; PDFs were
generated. Original /20mt/ and old outputs remain unchanged.

THREE scientific jobs remain LIVE; collect actual handles before source changes:
- 88009: outputs/material_tube_independent_20mt_viscous_r4_z16_t8000,
  same full array pulse with 0.4375 ns steps (time control).
- 85554: outputs/material_tube_independent_20mt_viscous_r6_z24_t8000,
  same pulse with 6x24 plasma cells, 6 vacuum cells / 12 modes, 0.4375 ns steps.
- 46579: outputs/material_tube_bias_only_viscous_r4_z16_t4000,
  matching bias-only control, 0.875 ns steps and declared C1=0.5, C2=2.
Last saved observations: 88009 at 600/8000 steps (378.14 s), 85554 at
200/8000 (302.81 s), and 46579 at 150/4000 (90.13 s); all nonterminal.
These are snapshots, not completion claims. Actual handles were re-polled live.
All preserve exact exterior vacuum and every original physical term, save every
accepted state, and use the same frozen source above. Do not restart from an
observation timeout or edit their source while live. The first two are tracked
by frc-vacuum.3.3; bias comparison remains in frc-vacuum.3.2. The completed numerical
failure/trajectory issue frc-vacuum.3.1 is closed; the full goal and parent remain open.

NEXT collect the bias-only run and compare its complete trajectory with the array
using examples.material_experiment_compare. Collect both refinements; use existing
actual material-basis sampling to measure B/J/velocity, geometry, heat and linkage
changes at common times/labels. The older material_coil_plot reconstruction expects
its older flat configuration, so do not feed these experiment configs to it unchanged.
No new verification framework or unit tests. Preserve the user-modified poster.

## 10ns animation and completed long-array trajectory (2026-10-16)

User requested denser web animation. The original completed two-loop20mT run
now has351display frames over0..3.5us at nearest10ns targets. Its actual solver
step remains1.75ns, giving8.75/10.5ns frame gaps; timestamps and that distinction
are visible. No temporal interpolation or scientific-state change. Read
examples/material_experiment_display.md and visualizations/material_fields/viewer.md.

The original URL127.0.0.1:8790/20mt/ now serves
outputs/material_experiment_20mt_10ns_display through the existing90454 server.
The same new view is at/20mt-10ns/; the previous12-frame view stays at/20mt-coarse/.
Dense local frame files keep the HTML around201KB, with only five complete frames
cached in the browser. All path vertices/current samples remain, with256 rendered
color levels for faster drawing. Previous/next, playback speed, play/pause, rapid
scrubs, signed components, color scales and final3500ns state were exercised in
the browser with no console errors/warnings. Native screenshotbrowser_view.jpg
was inspected. Full-page browser capture briefly produced a bad layout snapshot;
the saved screenshot was replaced by the correct native capture. No unit tests.

The actual long-array viscous run13616 also completed successfully: all4000steps,
3.5us,2456.917s. Source05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75;
outputs/material_tube_independent_20mt_viscous_r4_z16_t4000. Final middle/max radii
61.600/64.723mm, particle RMS25.318mm, Joule238.046mJ, numerical ion heat15.634mJ,
meanTe/Ti13.739/1.064eV. Max observed external-work residual/E0=4.5985e-5. This
completes a numerical trajectory, not convergence, physical closure or formation.
Final max drift/thermal ratio67.87 remains an important validity limitation.

ONLY scientific run84694 remains LIVE: the matched two-loop viscous control,
outputs/material_tube_two_loop_20mt_viscous_r4_z16_t4000;2700/4000 at last saved
observation. Collect its actual terminal result before source changes. Then
compare these matched runs and build the new experiment display if physically
interpretable. The current requested10ns animation deliberately retains the older
completed case; it does not silently substitute the new numerical treatment.
Keepfrc-vacuum.3.1/.3.2 and the full goal open. Preserve the user poster; no remote.

## Active paired pulse runs (2026-10-16)

Checkpointb72ae0bb7 saves the nanosecond timing label, complete shock comparisons,
optional numerical compression treatment and prior failed-run evidence. The new
long-array pulse has passed the previous1.101625us mesh-fold time: at recorded
1.40us (1600/4000steps,979.61s) middle/max radii34.476/34.943mm, Joule heat132.433mJ,
numerical ion heat10.900mJ, work residual/E0=5.2092e-6. This is an unfinished
trajectory, not a completed pulse or convergence/formation claim.

Two actual scientific jobs remain LIVE and must be collected before source edits:
-13616: outputs/material_tube_independent_20mt_viscous_r4_z16_t4000; nine-loop
 formation array plus independent bias;1600/4000 at the last saved observation.
-84694: outputs/material_tube_two_loop_20mt_viscous_r4_z16_t4000; matching two-loop
 control; initial accepted step collected. Both request3.5us with0.875ns steps,
 C1=0.5,C2=2 numerical compression pressure, full original physical terms and
 exact exterior vacuum. Both source05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75.

The matched control is necessary: comparing the new viscous half-dt array run
only with the earlier two-loop nonviscous1.75ns run would mix geometry and
numerical changes. Read MATERIAL_COMPRESSION.md. If both complete credibly, use
the existing physical comparison on these matched runs, include the separately
observed numerical heat, then build/inspect the actual3-D display. If either
fails, preserve its terminal evidence and address that concrete failure. Do not
restart live jobs from observation timeouts or substitute a partial prefix for a
completed experiment. No new framework or unit tests; keep the overall goal open.

The analytic prescribed-source comparison is in
outputs/material_source_geometry_comparison (PNG inspected, PDF generated).
Both centre fields are-19.9mT including bias, but at the initial plasma ends the
source fields are-0.638mT (two loops) and-17.142mT (long array). That is3.21% versus
86.14% of the centre magnitude. The same central Helmholtz bias is only3.207uT
at those ends. These are vacuum input profiles, not actual screened fields or an
equilibrium calculation. Hardware remains abstract prescribed waveforms.

The existing completed viewer stays at127.0.0.1:8790/20mt/ (server90454) with its
1.750ns integration label; its scientific data is unchanged. Preserve the user
poster. No remote/publication. The app goal record still reports paused; this
session is following the user's explicit continuation without marking it complete.

## Nanosecond integration and the strong-pulse compression limit (2026-10-15)

Research resumed after the usage refresh. The displayed completed two-loop20mT
trajectory uses2000 RK4 steps over3.5us: **1.75ns per solver step**, not250ns.
The12 selected display states are unevenly spaced by roughly0.25/0.5us. The
viewer now labels the integration step and display gap separately; its saved
scientific states are unchanged. The refreshed localhost8790/20mt/ view and timing
screenshot were inspected, with no browser warnings/errors. Server90454 serves it.

Both previous overnight jobs are terminal: the strong split pulse also failed
(613 accepted, attempted614, nonpositive/nonfinite species), while the higher-order
compact-vector trajectory completed2000steps. See MATERIAL_VECTOR_COUPLING.md for
its actual energy/virial and resolution differences; it is not spatially converged.

The full strong long-array pulse rerun with **0.875ns** steps reaches1.101625us
(1259 accepted, attempted1260) and rejects a folded material cell. This moves beyond
the original1.064us field/entropy explosion but does not complete the3.5us pulse.
A short0.4375ns continuation of the earlier prefix also reaches an end-cell fold;
it is diagnostic evidence, not an initial-to-end convergence run. Saved actual RK
stages, restricted EM spectra and folded geometry are in
outputs/material_independent_pulse_failure_analysis. The restricted spectrum does
not establish a simple electron-gyro timestep limit or full-system stability.

An optional numerical bulk pressure and matched ion heating are implemented in
material_bulk_viscosity and the experiment builder, defaultOFF. The current
complete wall-shock comparison11487 finished all50/100/200-cell trajectories,
with density/velocity/pressure L1 errors falling to0.962%/0.484%/2.956%. Read
MATERIAL_COMPRESSION.md and outputs/material_wall_shock_viscosity_b (inspected
PNG; PDF also generated). The finest energy/work residual is0.0523%. The complete
100-cell time control31436, outputs/material_wall_shock_viscosity_t4, cuts dt25ns
to6.25ns and lowers energy/work error0.105% ->0.000308%; physical shock errors
remain spatially limited. The physical20mT array
rerun is LIVE13616 in outputs/material_tube_independent_20mt_viscous_r4_z16_t4000,
source05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75;
0.875ns,4000requestedsteps, explicit numericalC1=0.5,C2=2. Collect actual terminal
results; do not restart on observation timeout or edit live sources. All previous
physical terms and exact exterior vacuum remain; no floors or mesh repair.
The first rejected one-radial-cell setup is preserved separately, with no evolved
state. No unit tests or new verification framework. The FRC goal andfrc-vacuum.3.1
remain open; physical closures, preparation, accuracy and formation are unfinished.

## Wrap-up: independent bias and formation coils (2026-10-14)

User requested wrap-up to conserve usage credits. No further scientific work is
being initiated. Remaining concrete work is filed as frc-vacuum.3.1 (strong pulse
failure/completion) and frc-vacuum.3.2 (complete physical comparison and display).
The overarching FRC goal remains incomplete.

The solver now supports independent prescribed coil groups, with actual summed
source fields/electric drive and separate linkage/work. The new configuration
keeps the100uT Helmholtz bias pair and adds nine100mm-radius loops spanning
z=-0.4..0.4m for a-20mT formation pulse. Initial gas/plasma and timing are unchanged
from the earlier two-loop20mT case; exterior matter/conductivity stays exactly zero.

The original RK4 long-array run49224 exited1 at attempted step609 after608 accepted
steps (1.064us),392.127s. exp overflow/nonpositive or nonfinite species stopped it;
all accepted states and frozen source69ef0ed6ad2e144ece369be2850d22940152fc0a02ebf7cb9b87f48783eb8a35
remain in outputs/material_tube_independent_20mt_r4_z16_t2000. At recorded1.05us,
middle/max radius34.663/34.980mm and heat60.961mJ; the subsequent final accepted
state already has explosively large field/current. These are numerical-failure
evidence, not credible predictions or completed formation. No floors/clipping.

The same physical case is now using the existing opt-in split method. Direct
external work remains observed; ambiguous instantaneous split response-work
integrals are null. This run is unfinished; weak-pulse success does not establish
strong-case stability. LIVE scientific jobs remain running:
- 70492: outputs/material_tube_independent_20mt_split_r4_z16_t2000, 200/2000 at last observation; terminal=False. Source 570bfe80c034da277c9aa5244aba3425f31bd66b7999b88f1bdd38e6c1668f2d.
- 89864: outputs/material_vector_free_r6_p5_q10_t2000, 1800/2000 at last observation; terminal=False. Source 5744fb58115ea903ba9a683959e4e47b30b331d75981d11bb672ae7e4fa743e6.
Poll these exact handles before changing any source. Preserve frozen closures
and every accepted state; never restart from a timeout. Last observed counts are
snapshots, not completion claims. Local execution can continue after this wrap-up.

The new display supports separate coil colors and centre-field labels. Explicit
--allow-partial shows terminal failed prefixes and their stopping reason. Builder
38088 exited0 for outputs/material_experiment_independent_20mt_partial_display;
3D PNG inspected. Browser controls and section PNG/PDF review are PENDING, and it
has not replaced the completed viewer at127.0.0.1:8790/20mt/ (server89915).
The comparison helper now supports group labels, particle RMS radius and independent
exterior-coil versus volume-current linkage. Its new path still needs execution
on complete runs; do not describe it as already verified. No unit tests were run.

NEXT collect70492 and89864, address the actual long-array failure if needed, then
use the existing comparison/display on complete physically interpretable runs.
Do not expand verification frameworks. Physical transport, ionization/preparation,
spatial accuracy and experimental calibration remain open. Research search found
additional primary papers but none was adopted as a new matched baseline here.
Preserve user-modified outputs/frc_showcase/poster.png; no remote/publication.

---

# Continue the physical FRC simulation

Read newest GOAL_STATUS.md and MATERIAL_EXPERIMENT.md. User priority: simulation
first, complete end-to-end and independent mathematical/experimental comparisons.
No unit tests, recursive verification frameworks or subagent delegation.

Completed bias-only,1mT and20mT pulses now distinguish magnetic response,
transient central compression, bulk extent and outer expansion. The20mT pulse
shrinks middle35 ->31.3557mm, then rebounds; final max radius68.6460mm. Particle
RMS radius25.6110mm vs27.4959mm control (both expand from18.7082mm). Heat85.5213mJ.
Strong cold-edge drift/thermal ratio116.93 and local accuracy limit physical claims.
Newest3D view localhost8790/20mt/ (server89915, Tab8 marked deliverable) includes
middle/max radius, actual full plasma/vacuum B and all J components. PNG/PDF,
comparison and browser screenshot inspected; controls checked without errors.

Completed corrected split runs45710/14252:1000/2000steps, same source
 d56d01a8aa762d4b58f241fb8abd3b5636278de928764c24497fdca1de1a0eeb.
Time-halving final J difference.38724%; weak Jr/Jz3.7249%/4.0554%, versus older
15.28%/14.99%. Fine work residual/E0 2.65e-10, entropy minima within2.91e-13.
This improves time integration but leaves spatial and strong-case validity open.
New step remains opt-in. Four failed attempts and the concrete bound-stagnation
fix are preserved. Existing physical plots perform comparisons; no new framework.

ONLY LIVE scientific job:89864 outputs/material_vector_free_r6_p5_q10_t2000,
older high-order full-vector trajectory. Poll exact handle; preserve its common
science source closure and all accepted states; never restart from timeout.
All new coil, control, strong-pulse, split and plot jobs are terminal.

NEXT prioritize physical simulations with independently prescribed bias and
formation coils, and defensible transport/ionization/preparation. Apply the new
integration to concrete physical cases before broad adoption; its weak-pulse
success does not validate nonconvex pressure branches or a laboratory experiment.
Collect89864 when it exits. Keep full goal and frc-vacuum.3 open. No remote or
publication; preserve user-modified outputs/frc_showcase/poster.png.

---

# Active continuation: full current coupling on moving material

## Current boundary: configured metre-scale biased pulse and 3-D view (2026-10-14)

Previous goal turn45221c28 was PROGRESS: actual moving vacuum and finite-coil
coupling, with two completed weak pulses. This continuation adds concrete
experiment inputs and completes a metre-scale opposing-pulse trajectory.
Read MATERIAL_EXPERIMENT.md and the material_experiment module companions.

`configs/material_tube_bias_pulse.json` describes a1m x15cm tube, Ar inventory
from0.1mTorr/300K, explicitly ASSUMED fully singly ionized/compact preparation,
R35mm,Lhalf350mm,nc8.96143e19/m³,centralTe/Ti10/1eV,eta1e-4ohm-m. Two prescribed
loops at r100mm,z±50mm supply100uT DC bias and a-1mT pulse. Rise/hold/fall1/.5/1us
starts at.25us; all2000steps/3.5us complete (20691 exit0). Actual full moving
current/fluid/species/vacuum terms remain; tube contact rejects at every stage.
No wall, neutral, ionization or preparation physics is silently supplied.

Peak Jtheta84.73kA/m², final Joule heat195.60uJ, source work156.49uJ. Plasma radius
ends54.35mm inside75mm tube; meanTe/Ti3.49554/.375728eV. After return to+100uT bias,
actual sampled plasma Bz spans-142.56 to+171.44uT. This is an exploratory response,
not converged FRC formation or confinement. Maximum electron drift/thermal ratio
5.95 and final cold-edge ion-ion Coulomb log.0159 expose physical validity limits.
Maximum external/source-work residuals3.88e-7/4.12e-7 do not erase those limits.

The completed12-state3-D view is outputs/material_experiment_display/index.html,
served at http://127.0.0.1:8790/ (server89915), browserTab8 marked deliverable.
It includes actual exterior B, prescribed coils, tube, moving plasma and all J
components. Scientific3-D PNG/PDF, signed Bz/Jtheta/flux sections and browser_view.jpg
were inspected. Playback, all current components, zero-current initial frame,
shared/per-state color scales, rotation, zoom and reset exercised; browser logs
show no warnings/errors. All scientific states remain unchanged. The old compact
view at8789/server55641 and the user poster remain preserved.

The6x12 weak coil refinement6211 also exited0. Common material-label field/current
comparisons against4x8 show J differences30.26%,8.56%,60.38% at quarter/half/final
duration. Final azimuthal J alone differs59.43%; total B difference23.62% after
the source is off. Final heating differs0.98%; integrated agreement does not
establish local accuracy. Independent coil-linkage discrepancy.3791% ->.3074%.
Read outputs/material_coil_spatial_comparison_q10. The common10-point quadrature
samples actual basis fields; no quadrature-cloud interpolation or state fitting.
The preliminary6-point result is preserved. No unit tests or new test framework.

LIVE scientific handles, authoritatively polled at this boundary:
- 63566: outputs/material_tube_bias_only_r4_z16_t2000. Same metre-scale assumed
  state, no pulse.525/2000,307.1s,RSS208MB at observation.
- 42299: outputs/material_coil_free_r8_t2000. Next actual weak-pulse spatial grid;
  775/2000,884.6s,RSS287MB at observation.
- 89864: outputs/material_vector_free_r6_p5_q10_t2000. Older compact high-order
  full-vector trajectory;925/2000,5250.2s,RSS492MB at observation.
Poll exact handles; do not restart from timeouts or edit their source closures.
Terminal display38726 and spatial figure44644/13723 exits0 collected.

NEXT compare the complete metre-scale pulse with its bias-only control using
actual fields, current, heat and motion; collect the third weak-pulse spatial
trend and older high-order trajectory. Address numerical local-current limits
and physical transport/edge/ionization, then compare meaningful biased/pulsed
conditions and experimental observables. MATERIAL_EXPERIMENT cites the primary
HFRC-F formation measurements: main-phase ionization/preparation matters, so the
current fully pre-ionized calculation is NOT an experimental calibration.
Keep frc-vacuum.3, transport, formation and the full goal open. No remote/publication.

## Current boundary: finite coils drive the freely moving plasma (2026-10-14)

This continuation of836bcad9 advances the actual simulation. MaterialMagneticVacuum
follows the freely moving side and both ends with zero exterior plasma matter/
conductivity; its condensed magnetic energy, metric rate and shape force enter
MaterialVacuumCurrent. MaterialVacuumSpace matches the same poloidal flux trace
throughout the interface. MaterialCoilCurrent supplies analytic finite-loop
fields, stage-time electric drive, motional terms, Lorentz force/torque and two
independently accumulated reciprocal source-work budgets. Full physical masses,
Hall/current inertia and stress, rotation, species entropy/pressure and heat
remain. Read MATERIAL_VACUUM_COILS.md and module companions before changes.

Three complete nonaffine moving-vacuum trajectories compare against analytic
physical dipole fields, Maxwell traction, energy and power. Joint refinement4/8/16
reduces maximum field error9.06% ->1.38% ->0.266%, traction4.60% ->0.542% ->0.0806%,
and power3.55% ->0.347% ->0.0134%. An independent physical sampler also converges.
This is a mathematical baseline, not experimental calibration.

The complete unforced penetrating-flux trajectory40493 exited0:2000steps/7.12us,
energy change3.00e-6, all-space virial error3.43e-7. Initial B/J errors19.75%/14.31%
remain coarse. Its frozen archive retains the temporary module name
material_vacuum_expansion; that ORIGINAL nonmagnetic driver was restored unchanged.
The new permanent entry is examples.material_penetrating_expansion. Both startup
failures are preserved with their frozen source and no physical trajectory.

Both complete10uT finite-coil pulses,1000/2000steps, exited0 (30464/26556).
They include rise, hold, fall and source-off evolution through7.121769us.
At1.780us, time halving changes total B/J by0.00193%/0.1013%; final changes are
0.624%/1.498%. Componentwise final J(r,theta,z) differences15.28%/0.150%/14.99%
and Btheta16.24% reveal unresolved weak late-time currents despite the dominant
component's agreement. Maximum external-work residual/initial energy improves
7.60e-9 ->3.54e-9. Final Joule heat5.49395nJ, source work5.53643nJ. During the
plateau, exterior flux at actual coil locations differs0.3791% from volume-current
linkage; neither is adjusted to fit the other. Inspected figure/data:
outputs/material_coil_comparison/material_coil_pulse.png and comparison.json.

The heated upwind time control56380 also exited0. Both full trajectories preserve
the initial entropy lower bound, but final B/J/drift differences8.74%/14.36%/19.89%
are not time converged. Maximum energy drift improves only2.21-fold. The matched
pressure force's upstream switch at zero current is a candidate source of weak
late-current sensitivity, not a proved sole cause. No entropy/heat/transport floor
or fitted correction was added. Read MATERIAL_THERMAL_COUPLING.md.

LIVE scientific handles (authoritatively polled after the above completions):
- 6211: outputs/material_coil_free_r6_t2000, same complete pulse on6x12,p3,q6,
  vacuum4/modes6,2000steps. At latest observation575/2000,361.8s,RSS208MB.
- 89864: outputs/material_vector_free_r6_p5_q10_t2000, older compact full-vector
  initialized-flux run; at525/2000,2970.6s,RSS492MB. Initial B/J0.432%/0.844%.
Poll exact handles for actual terminal exits; do not restart on observation timeout
or edit their frozen source closures. Keep all accepted states and failures.
The old3-D compact-vector display remains at127.0.0.1:8789 (server55641); it does
not display the new coil pulse. Other old scientific live notes below are history.

NEXT finish the6x12 coil pulse and compare actual fields/current, boundary motion,
heat and coil linkage with4x8; inspect both dominant and weak current components.
Resolve the measured temporal/local-current limits before stronger-waveform
interpretation. Then advance physical transport/ionization and useful prescribed
bias/formation studies, with experimental comparisons where applicable. No new
verification framework or unit tests. The present centimetre-scale fully ionized
fixture is not the user's proposed metre-scale neutral-fill experiment. Spatial
and temporal convergence, heated-edge/transport validity, formation, experimental
controls, stability and confinement remain distinct unfinished requirements.
Keep frc-vacuum.3, its parent and the full FRC goal open. No remote/publication.

## Current boundary: full-vector moving plasma and updated 3-D display (2026-10-14)

The preceding checkpoint5e61a3eb was progress: actual species heat and entropy
transport. This continuation evolves all three axisymmetric magnetic/current/
flow components, poloidal flux and actual material rotation with finite current
inertia, Hall, current stress, both species and heat. Read
MATERIAL_VECTOR_COUPLING.md. A complete4x8,p3,q6,2000-step trajectory finishes
7.121769us. Its torque generates8.30m/s peak toroidal motion. Energy drift2.33e-6,
independent virial error2.01e-7, axial/angular momenta remain near roundoff.
Initial B/J errors19.6%/35.9% are large: this is NOT locally converged.

The new3-D display in outputs/material_vector_display/index.html uses five actual
saved states, interpolated full-B paths, all current components and the moving
boundary. It is served at http://127.0.0.1:8789/ (server55641) and opened in browser
Tab7. Time/component/scales/playback/rotation/zoom/reset were exercised; no browser
warnings/errors. PNG/PDF and browser_view.jpg inspected. Initialized flux decays;
this does not demonstrate FRC formation, confinement or3-D stability. Old viewer
and user poster are preserved.

The first complete upwind6x12,p5,q10 thermal run60494 exited0, preserving the
initial entropy lower bound throughout7.121769us; max energy drift1.12e-6 and
virial error1.48e-7. First-order numerical diffusion remains explicit. Physical
edge validity remains unresolved (minimum Coulomb log1.026). Earlier failed skew
run26494/stage187 remains preserved. No entropy/heat/log floor or collision
removal was introduced. A time-halved full control is still live below.

The nonresistive8x16,p5,q10 current run43078 also exited0: initial B/J errors
0.0683%/0.2096%, max energy drift1.01e-8, virial error5.94e-11. Nevertheless final
magnetic/current energies31.80/53.94nJ differ substantially from degree3's
15.41/63.98nJ; peak J changes14.28 to6.07kA/m². Better initial accuracy and good
conservation do not establish full spatial convergence. Read the inspected
outputs/material_current_p5_comparison figure and comparison.json.

LIVE scientific handles:
- 56380: outputs/material_thermal_free_r6_p5_q10_t2000_upwind, time-halved control.
- 89864: outputs/material_vector_free_r6_p5_q10_t2000, actual full-vector6x12,p5,q10.
  Independent INITIAL B/J errors0.4320%/0.8437%; acceleration RMS0.0018903 of C/R.
  This larger coupled run takes substantially longer than the coarse trajectory.

Re-poll exact handles, collect actual terminal exits and inspect summaries; do not
restart on observation timeout or edit their source closures. Both preserve every
accepted state and frozen source. Scalar-gradient toroidal rows are exactly zero;
MaterialDualEntropy now admits them without changing meridional transport.
Terminal scientific20443,43078,60494 exits0 collected; display90906/6534 and
comparison84237 also terminal0. Other old logs below are historical.

NEXT finish the two live trajectories and compare actual thermal time changes
and full-vector fields/motion. Upwind energy drift needs its own measured time
assessment; do not borrow the old skew method's fourth-order result. Then advance
penetrating moving vacuum/finite coils and physical transport/ionization. The
current vector boundary enforces B=0 on the free surface with exactly empty
exterior. It is not the finite-coil formation model. Preserve the older fuller
fixed-footprint coil solver. Keep frc-vacuum.3, parent and full FRC goal open.
No unit tests or new verification framework. No configured remote/publication.

## Current boundary: finite-mass current coupled to the free material (2026-10-14)

This advances c4c32463 from ideal magnetic forces to independent current/field
coordinates evolving with actual material motion. MaterialCarrierCurrent retains
physical electron/ion masses, Hall/self-drift, moving-density inertia, vorticity,
magnetic force and drift stress. Read MATERIAL_CURRENT_COUPLING.md and the
visually inspected outputs/material_current_comparison/material_current_expansion.png.
The invariant sector is toroidal B with meridional flow/current, equal constant
species entropies and eta=0. Both ends and side move freely; exterior matter,
conductivity and toroidal B are exactly zero. It is a mathematical component,
not a replacement for the fuller fixed-footprint finite-coil solver.

Two complete prescribed-expansion trajectories follow the independent finite-
inertia law (physical B amplitude 1/6 at doubled scale, versus ideal 1/4).
Final B/current errors improve 4.88/2.70-fold with spatial refinement. Five
free-run attempts include four completed 7.12177us trajectories and one preserved
failure after553 accepted steps. Original fine Faraday reconstruction exceeded
its unchanged1e-10 limit. Refining the original mixed linear rows resolves that
failure; no state/energy correction or tolerance relaxation was introduced.

Independent initial data exposed magnetic-only projection current-energy errors
+84.7%/+29.2%. A physical magnetic-plus-carrier error projection reduces those to
+3.49%/-0.152%, but increases initial B error. The finest initial J/drift errors
remain16.2%/25.4%. New free trajectories conserve energy within6.44e-10 and satisfy
the independent virial within3.56e-12(coarse)/3.22e-11(fine). Yet final magnetic
energy changes45%, current energy22%, and peak current4.5-fold between4x8/8x16:
THE FREE TRAJECTORY IS NOT SPATIALLY CONVERGED. Global conservation does not prove
accurate local fields, pointwise current continuity or physical edge transport.

Conditional Coulomb samples flag eta=0 (34–190 collision times over the run) and
large edge drift (final drift/electron thermal speed up to24). These conditional
Maxwellian estimates do not validate transport at large drift. No laboratory
prediction, FRC formation, stability or confinement is established. General
species entropy/heating, poloidal fields and moving vacuum/coil coupling remain.

All scientific sessions15237,1742,50963 and batch62098 have terminal exits0
collected; original failure94090 exit1 and observers95334,43866 exit0 collected.
No worker is live. Do not restart them. Frozen sources and every accepted state
remain in outputs/material_current_*. Figure PNG inspected; PDF/comparison JSON
saved. No unit tests or new verification framework. User poster remains untouched.

NEXT improve the actual field/current spatial representation using existing
polynomial degree/quadrature controls (currentlydegree3/q6; candidate5/q8), then
compare a complete free trajectory and necessary time control. Do not replace
simulation work with another reporting framework or claim convergence from
energy alone. Then extend general species entropy/heating and poloidal/moving-
vacuum/finite-coil physics. Keep frc-vacuum.3, parent and full goal open. The old
interactive3D viewer remains unchanged. Earlier boundary blocks are historical.

## Current boundary: magnetic forces on the moving material (2026-10-14)

The previous goal turn made progress in05b4abd4 by evolving the material fluid.
This turn adds actual full-vector ideal magnetic forces and angular torque,
including field winding by material rotation and the full geometry derivative.
Read MATERIAL_MAGNETIC_COUPLING.md and the inspected figure in
outputs/material_magnetic_comparison/. This is component progress, not completed
finite-mass moving plasma or a replacement for the full finite-coil solver.

The reference field vanishes at the material boundary; the exterior is exactly
empty with zero B. General nodal motion is evolved, with no prescribed shape,
force fit, diffusion or energy repair. This ideal fixture has initialized closed
flux, no coils, and no resistive formation. Current inertia/Hall/species-relative
transport/Joule/exchange are absent from this moving path and remain mandatory.

Four complete7.12177us trajectories (initial r4, then4/8/16 with corrected current
observer) finish; actual scientific95584 and batch18833 exits0 were collected.
Plot99255 exited0; PNG inspected. No worker remains live. Do not restart them.
Current reconstruction was corrected to include internal field-jump contributions
through weak Ampere; cellwise curl alone has a nonvanishing refinement error.
The projection is a diagnostic, not a finite-mass charge/current discretization.

Independent initial force errors0.01695/0.004332/0.001064 of C/R; weak-current
errors0.1382/0.04542/0.01163 of Bp/(mu0 R). B errors improve only about2x per
mesh doubling, not4x. Complete final position/velocity/current differences
improve about4x; B differences2x. Native energy change<2.89e-14 and independent
virial inertia error<6.38e-10; initial magnetic-energy geometry error0.174% remains
at16/32. Actual torque generates4.54m/s peak toroidal motion. Cold dilute edge
and projected J/n show strong validity/reconstruction limits; no lab claim.

NEXT implement independent magnetic/current coordinates on the moving geometry,
with full finite-mass energy/stress and Hall/species transport, then true finite
coil and vacuum flux traces/source work. Do not make further ideal-only studies
or reporting frameworks the next milestone. The current reconstruction reveals
why raw cellwise curl/J/n is not a valid coupling. Keep frc-vacuum.3, parent and
full goal open. Preserve the existing full finite-coil solver, all historical
sources, and the user's modified poster. The older interactive3D viewer is
unchanged; formation/transport/UI/experimental acceptance remains unfinished.
Earlier sections below are historical and superseded.

## Current boundary: moving carrier/thermal fluid component (2026-10-14)

The previous goal turn made progress in95bd276: complete finite-coil/sealed-end
coupling and seven full pulses. This turn adds a necessary moving fluid
component, not a replacement or completed extension of that electromagnetic
solver. Read MATERIAL_CARRIER_EDGE_RESULTS.md and the inspected
outputs/material_vacuum_comparison/material_vacuum_expansion.png.

MaterialTwoSpecies carries fixed quadrature particle weights, two physical
entropies, canonical meridional/angular momentum and material rotation angle.
Density uses actual changing cylindrical volume. Compact reference density and
temperatures vanish together at the moving surface; no exterior fluid/floors.
Consistent physical mi+me inertia, pressure virtual work and position-dependent
angular inertia/centrifugal force evolve general Q1 node motion. Outer motion
is not prescribed or constrained to affine shape. The older P0 MovingPlasma
solver is not reused as dynamics. No unit tests, new framework, diffusion,
energy repair or acceptance relaxation was introduced.

Three spherical and two rotating elongated7.12177us trajectories complete;
a sixth complete rotating trajectory stores/compares material phase. Actual
sessions90694,87440,12290 and figure observer95156 exited0, all collected.
Initial directr4 setup exited1 before stepping on a name collision, preserved.
No worker remains live; do not restart these runs. All accepted states/source
closures remain in outputs/material_vacuum_*. The final angle-aware example
uses the current source; earlier five frozen sources omit the passive angle.

Independent physical-vacuum Euler reference: closed-form spherical expansion
and separately integrated rotating affine axes. Roughly4x error reduction per
mesh doubling. Finest final surface errors40.0um(sphere)/24.5um(rotating),
velocity errors0.0348%/0.0308% of stated speed scale, density errors~0.083%.
Native energy drift <4e-14; initial polygonal mass remains0.2408% low at16/32,
so conservation is not physical accuracy. Final phase error0.000100835rad.
B=J=0 and equal species temperatures make current/Hall/Joule/exchange terms
vanish in this explicitly limited fixture. No heated edge, ionization or
laboratory-validity claim. Rotation/volume are real evolved solver states.

frc-vacuum.3 and full goal remain open. NEXT connect full magnetic/current
physics to the actual moving material state, rather than add more zero-current
or sphere/refinement campaigns. Include changing geometry, toroidal winding,
physical current kinetic energy, Hall/advection/stress, species relative flux,
finite-coil mutual/source work and a magnetic-only exterior. The new
axisymmetric deformation helper is preparatory, not verified field coupling.
A carrier-weighted drift space is a candidate to avoid artificial J/n floors,
not an adopted closure. Read SURFACE_COMPATIBILITY before edge heating; the
zero-current physical-vacuum profiles do not solve current-carrying thermal
compatibility or kinetic/sheath/ionization validity. Retain full finite-coil
solver, all archived evidence and user poster. Full formation/UI/experimental
validation, long-contact accuracy and remaining transport gates stay open.
Earlier boundary blocks below are historical and superseded.

## Current boundary: finite coils and sealed finite ends (2026-10-14)

This goal turn made physical simulation progress beyond23b5d08. The full
finite-mass moving plasma now couples to actual circular coils through a
magnetic-only exterior surrounding the side and both ends of a finite tube.
No periodic axis or periodic coil array; exterior has exactly no matter or
conductivity. Full current inertia, Hall, transport/stress, species, heat,
mutual magnetic energy and prescribed-source linked-flux work remain.
No unit tests, new verification framework, floors, numerical diffusion, energy
repair or residual relaxation. Hardware remains prescribed waveforms.

Read FINITE_CYLINDER_RESULTS.md and the visually inspected
outputs/finite_coil_comparison/finite_coil_response.png. Five full acoustic
periods, four static dipole comparisons and seven full24us coil pulses finish.
Actual scientific sessions10140,15748,37631,59968,48323,98433,71975,74225
are terminal with exits0 collected (first vacuum case direct execution0).
Initial coil observer setup9947 exited1 before stepping and is preserved;
figure observer69616 exited0. No worker remains live. Do not restart them.

Weak-conductivity independent current error falls3.35%→0.183% with field
refinement; time halving improves energy/work error16-fold. Full10uT Spitzer
runs on8/16,16/32,24/48 fields retain actual motion and both temperatures.
The last refinement changes peak current0.3314% and heat0.000993%; mean-field
change is4.34e-7 of imposed center peak. Finest peak J2081.11A/m2,heat0.554688uJ,
peak speed0.02019m/s. At6us mean plasma B is28.79% of coil-only B averaged
over the same volume. The100uT case crosses weak-field/drift review thresholds;
10uT cases stay below them, which is not complete transport validation.
Fixed exterior12-cell/8-mode linkage discrepancy remains~0.075%, not jointly
converged by plasma refinement alone. Sealed end caps omit open-end/sheath loss.

frc-vacuum.2 is complete for the bounded finite-coil/sealed-end model and its
measured comparisons. The full goal and parent frc-vacuum stay open. NEXT:
frc-vacuum.3, a physical moving plasma edge retaining this complete current/
species model. Read SURFACE_COMPATIBILITY.md and MOVING_PLASMA.md first. A
smooth continuous-B,p=0 edge with finite eta J² has incompatible thermal traces;
select and justify an actual carrier/thermal edge model. Do not remove edge
heat, taper eta, add an atmosphere/floors, or substitute an older incomplete
solver. Implement an actual moving-interface trajectory against an independent
expansion/traction baseline, then finite-coil compression. No new framework.
Transport/ionization, strong-field applicability, long-contact accuracy and
converged formation remain. No experimental calibration,FRC formation,3D
stability or confinement is established. The older3D viewer remains its old
contact snapshot. Preserve the user's modified outputs/frc_showcase/poster.png.
Earlier boundary blocks below are historical and superseded.

## Current boundary: magnetic vacuum coupling implemented (2026-10-14)

The previous goal turn made progress: commit 49b0102 completed the driven-wall
pulse comparisons. This turn implemented an evolving plasma/vacuum magnetic
boundary and completed five full 24 µs runs. Read VACUUM_CYLINDER_RESULTS.md and
the inspected outputs/vacuum_pulse_comparison/vacuum_coupled_response.png.
All scientific sessions 66738,10198,15544,13948,44808 and plot observer 67929
are terminal with actual exits 0 collected. Do not restart these trajectories.

VacuumCylinderGeometry/VacuumPlasma/VacuumInduction/VacuumGauss retain full
moving finite-mass current/species dynamics. Outer flux is free. The exterior
contains no plasma density/conductivity; radial K-Bessel elimination uses the
same polynomial axial trace and includes induced energy to infinity. An ideal
cylindrical source sheet drives the plasma, with linked-flux work and the nonzero
source/induced cross energy included. No new framework, unit tests, numerical
diffusion, artificial atmosphere, state repair or tolerance relaxation.

At the 100 µT source peak, mean plasma B is 51.914 µT at fixed resistance and
57.625 µT with evolving resistance. The latter yields 8.7839% more heat in this
new prescribed-source configuration. Joint 8/16→16/32 spatial refinement gives
current/vacuum space-time errors 0.001085%/0.0003213%, improving 5.50/29.60-fold.
Time halving improves energy error 15.97-fold. Global errors do not describe
startup uniformly: the zero-field coarse control has 4.31% current error at its
first 150 ns endpoint. The axial seed is a mathematical fixture, not formation.

frc-vacuum.1 is complete for this fixed-footprint periodic coupling only.
NEXT physical work: frc-vacuum.2, finite coil fields through actual axial end
boundaries. Keep the actual new plasma equations and exact empty exterior;
reuse existing coil/open-vacuum machinery only with compatible trace/work maps.
Do not substitute a periodic coil array or an older incomplete plasma solver.
No finite-end implementation exists yet. Avoid another reporting milestone.

The complete FRC goal remains active. Moving free interfaces, actual coils/ends,
transport/ionization/validity, frc-7mm.7 long-contact accuracy and converged
formation remain. Preserve frozen evidence and the user's modified poster.
Earlier sections are historical; their live-worker notes are obsolete.

## Current boundary: driven component complete; vacuum coupling next (2026-10-14)

Previous goal turn bf56f1e made physical progress by evolving the missing wall
current and completing two pulses. This turn completed the full16-cell spatial
comparison and corrected Spitzer pulse. Read WALL_CURRENT_RESULTS.md and the
visually inspected outputs/wall_current_comparison/wall_current_refinement.png.
All scientific workers are terminal; actual sessions14381/65413 exited0, as did
plotting observer97286. Do not restart these completed trajectories.

Four full24µs runs now retain physical wall-current inertia and the same current
in all forces, heat and source work. Fixed eta on8→16 cells at320 steps reduces
space-time current error5.992-fold to0.003419%, peak Bz error6.831-fold to
4.10337e-5 of drive, and wall-power discrepancy15.329-fold to2.58659e-5W.
160→320 time steps reduce energy drift15.9919-fold. Global pulse norms are not
uniform accuracy: the first75ns current RMS error remains1.26% of its small
instantaneous reference; no pointwise current time-convergence claim is made.

The complete corrected Spitzer run increases eta38.3182→74.0081µΩm, gives
8.295722µJ thermal gain (17.38% below fixed eta), final maxJ0.73555 versus
52.11797A/m², and peak speed0.04929m/s. These are physical control differences,
not numerical errors. Weak-field/drift ratios pass review thresholds, not full
closure validity. No unit tests, new framework, numerical diffusion or repairs.
frc-goal.1 is complete for its bounded prescribed-flux component only.

NEXT physical solver work is frc-vacuum.1: replace fixed outer flux with an
evolved magnetic vacuum/source interface. Exterior gets no fluid density or
conductivity. Existing exterior machinery uses different traces, so couple
its actual field/work maps consistently; do not merely splice representations.
Free outer flux also makes the harmonic/increment cross energy nonzero; retain
it and source work. Begin the fixed-footprint coupled pulse, then continue the
parent's free interface, finite coils and physical end boundaries. This step
has not been implemented yet. Avoid another reporting/replay milestone.

Full FRC goal, transport/ionization/validity and frc-7mm.7 long-contact accuracy
remain open. Filled periodic rigid-cylinder results do not prove formation,
3D stability or particle confinement. Preserve frozen evidence and user poster.
Earlier live-worker notes below are historical and obsolete.

## Current simulation: prescribed flux and evolving resistance (2026-10-14)

Three complete 24 µs driven-plasma runs now finish: two fixed-resistance grids
and one evolving-Spitzer-resistance control. Read DRIVEN_FLUX_RESULTS.md and
outputs/driven_flux_comparison/driven_flux_response.png (visually inspected).
Actual sessions 4859,4781,70772 and plotting observer25683 all exited0; no
scientific worker remains live. Do not restart those calculations.

New AxialFluxDrive and DrivenModalGauss supply actual loop-voltage forcing,
source-time-dependent stages and external work. The pulse prescribes enclosed
flux pi R² Bbar, not wall Bz. R0.035m,L0.24m,n1e21/m³,Ar+,initial Te5/Ti1eV;
peak mean B0.1mT for12µs followed by12µs relaxation. Full momentum, current
inertia, Hall, stress, pressure, Joule heat and Coulomb exchange remain.
No numerical diffusion, atmosphere, floors, unit tests or new verification
framework. Saved accepted states/stages and frozen source closures remain.

Fixed eta3.83182e-5ohm m,160 steps,8→16 radial cells: maximum Bz RMS error
falls2.79186e-4→4.82744e-5 of drive amplitude; space-time current error remains
28.06%→21.14% of the continuum reference RMS. Heating is7.87%→4.47% low.
Independent finite-inertia Bessel diffusion retains nonzero wall current;
the old production current basis forces Jtheta(R)=0. Discrete-source versus
reconstructed wall-Poynting power discrepancy falls0.131842→0.0747768W.
Native energy closure cannot establish accurate current or boundary work.

The full evolving-resistance run cools electrons/heats ions to3.02134/2.97870eV,
raising eta to74.0081µohm m. Same-grid thermal gain is16.61% lower than fixed eta;
final maximum current0.73558 versus52.11798A/m². These are model/control
comparisons, not numerical errors. Peak Hall parameter0.01052 and ion drift/
thermal speed0.03093 pass existing review thresholds, not full validity proof.

NEXT physical work: correct the missing driven-wall current degree of freedom,
including its finite current inertia and paired work/heating. Do not add
compensating heat or alter current only in diagnostics. Use the completed
pulse as the mathematical baseline; avoid more reporting/replay infrastructure.
frc-goal.1 stays in progress until driven current/heating are credible.

These filled periodic rigid cylinders have no solved vacuum, finite coils or
real ends and do not demonstrate FRC formation. Full FRC goal, vacuum/transport
scope and frc-7mm.7 long-contact accuracy remain open. Prior hyperdiffusion
remains rejected. Preserve all archived evidence and the user's poster.
Earlier sections below are historical; their live-worker notes are obsolete.

## Current simulation boundary (2026-10-14)

Completed physical progress: three full100us moving temperature contacts and
four full12us collisional magnetic decays. No unit tests or new verification
framework. All scientific workers are terminal; actual exits were collected.

Nonuniform advection/Coulomb exchange matches independent characteristic
relaxation. Axial refinement reduces Te/Ti errors50/55.7-fold; time refinement
reduces energy drift15.99-fold to3.93322e-11. Read ADVECTED_COULOMB_RESULTS.md.
Sessions38573/39308/91750 completed; frc-transport.1 is closed for this component.

New SpitzerWallPlasma evaluates weak-field Z1 resistance from actual evolving
Te/n, paired with electron Joule heat and electron-ion exchange. Cooling raises
eta73%; magnetic decay follows an independently integrated variable-resistance
Bessel reference. Original8/80 and16/80 runs completed, sessions61936/53684.
Final B errors3.27565e-5/8.07015e-6 of remaining signal; current error does not
mesh-converge. Read SPITZER_DECAY_RESULTS.md and the inspected comparison PNG.

Original8/160 rejected step17 on strict floating wall-trace cancellation;
session58978 and targeted observational probe92257 are terminal with the
failure preserved. New opt-in LiftedModalGauss adds independent magnetic
coordinates before one lift, retaining all equations, scaling and bounds.
Its8/80 and8/160 full runs complete, sessions99314/35233: energy drift decreases
15.9858-fold to2.91367e-5 of initial field/current energy. Old/new80-step final
vector-B difference is1.90935e-10 of initial Bz. B/J errors do not time-converge;
finite-amplitude/spatial/numerical contributions remain. frc-transport.2 closed
for the bounded implemented component, not complete transport validity.

NEXT frc-goal.1: drive the collisional plasma with a prescribed axial-field
pulse. Couple induced-field Faraday forcing and external work explicitly;
merely changing the harmonic bias is insufficient. Begin with the filled
cylinder/fixed-resistance diffusion baseline, then evolve resistance, Joule
heat and motion. Hardware stays prescribed waveforms. Use existing tools and
actual end-to-end physics, not another replay/reporting framework.

The closure is conditional: weakly magnetized Z1 Maxwellians, leading-me/mi
Joule heat, no thermal forces/conduction/viscosity/neutrals/ionization. It is
not a complete magnetized collision model. Exact vacuum, physical coils/ends
and credible formation remain open. The hyperdiffusion candidate remains
rejected for long-time accuracy (HYPERDIFFUSIVE_CONTACT_RESULTS.md), and
frc-7mm.7 stays open. Preserve all frozen evidence and the user's poster.
Earlier sections are historical; their live-worker notes are obsolete.

## Current simulation run (2026-10-14)

Read COUPLED_FLUID_REFINEMENT_RESULTS.md. F32 audit and full cross-grid
comparison/replays are complete. Do not rerun them. The user specifically
prioritizes simulation and direct mathematical/experimental comparisons over
additional verification frameworks; AGENTS.md records that preference.

Actual full F32/64,B32/64 evolution runs in session78724, parent PID33604,
child33608, via the existing trajectory_control_campaign. Target100us/60 steps,
epsilon=-1e-4, original bounds. Poll that handle, collect its actual terminal
exit, then compare saved physical evolution with the independent reference.
Campaign: outputs/trajectory_fluid32_campaign; states:
outputs/trajectory_fluid32_t60_minus. One heavy worker at a time. No new physical
source has changed, and no further report/replay wrapper is needed for this run.
All full-goal gates remain open; older pending-comparison notes are historical.

## Latest boundary: requested signed-current / magnetic display (2026-10-12)

Use outputs/current_state_display_v3/index.html and current_state_3d.png.
Complete build/replay is byte-identical; PNG visually inspected,560 inputs
hash-checked. Signed current selection and meridional Jz/Bz-departure slices
add detail without changing the four physical states or magnetic paths.
Latest accepted step remains17 at28.3333us; the new diagnostics do not evolve
time. Mac locked, so browser interaction remains unchecked. Read
CURRENT_STATE_DISPLAY_V3.md. Optional localhost presentation server60135
uses127.0.0.1:8788; it is not a scientific worker.

NEXT scientific work: COUPLED_FLUID_REFINEMENT_PROGRESS.md records completed
transfer/F16/F32 probes and F16 physical audit, each fully replayed. Finish
the unlaunched F32 physical audit/replay, then the prepared but unexecuted
physical cross-grid comparison under a new guarded wrapper. Session30234
is expired; saved replay/terminal evidence passes, but no new actual exit
was observed. Do not restart it. No further physical job ran for the display.
All goal gates remain open; keep original bounds and use end-to-end checks.

## Latest continuation: completed entropy budget and scalar-space refinement

Read SPECIES_INTERFACE_VERIFICATION.md. All workers are terminal with actual
exits0 collected; 26 raw files and all scientific metrics/report PNG fully replay.
The whole weak entropy identity is now accounted for to1.78627e-12 relative.
Most of the saved strong-source discrepancy comes from unresolved fluid-space
representation, with a smaller mass-map effect. Holding current fixed, original
electron source error falls55.0%→8.97%→0.0292% on F8→F16→F32. This does not
demonstrate a coupled growth cure. All original physical sources/states/bounds
and the user poster are preserved; no unit tests or solver changes.

Final preservation record: outputs/species_budget_diagnostics/checkpoint_v2.json
(complete 200-source union; prior aggregate record retained).

NEXT COUPLED_FLUID_REFINEMENT_DESIGN.md specifies complete unchanged-RHS controls
after exact subdivision and independent verification of the full archived
directions, including dependent-axis transverse fields and both species. Keep
B/current fixed and the original normalization, use both signs/three amplitudes,
retain failed original checks, and complete full fresh replays. Then select full
time evolutions and local energy controls. All goal gates remain open; do not
restart old jobs or replace the solver mass with a diagnostic matrix.

## Requested display update (2026-10-12, after d4cff11)

Current visualization: outputs/current_state_display_v2/index.html and
current_state_3d.png. The PNG is inspected and full HTML/PNG/manifest generation
replays byte-identically. New latest-audit status, camera presets and keyboard
controls use the same four saved states ending at step 17, 28.3333 µs. No solver
was run or advanced. Browser interaction remains unchecked because the Mac was
locked. Read CURRENT_STATE_DISPLAY.md for evidence and display limits.
The full scientific next action remains SPECIES_INTERFACE_DIAGNOSTIC_DESIGN.md;
do not restart completed physical workers. All acceptance beads remain open.

## Latest continuation: physical spatial contributions and current flux (2026-10-12)

Progress from fa597fc: the complete original RHS is now decomposed and freshly
replayed on all12 saved signed tangent states. All96 separated induction solves
pass the original1e-10 bound. Full native fluid sums match bitwise; numerical
superposition errors in magnetic/current responses remain explicitly measured.
Independent physical q16/20 and continuum-term audits replay completely.
Read SPATIAL_PARTS_VERIFICATION.md and SPATIAL_PARTS_DESIGN.md.

Electron pressure drives the saved-error Btheta response (+1.61542e6/s physical
directional contribution), partly opposed by Hall (-6.64587e5/s). These are
signed component projections, not eigenvalues or an established cause. The
smooth independent mode exposes large cancellation among accurate individual
responses; its existing net-current error remains visible.

The fresh one-sided current audit reproduces all4 raw records and verifies
cellwise physical Gauss identities to5.3891e-16 relative to absolute face flux.
Saved-error R||divJ||/||J||=0.50247 and normalized radial jumps=0.05605, versus
0.00037794 and2.06491e-8 for the independent mode. Normal wall current is near
roundoff. Global cancellation does not establish local charge consistency.
Strong interior entropy/density discrepancies are about55%; omitted interface
terms prevent a causal conclusion. No term or state has been repaired.

All scientific jobs are terminal, actual exits0 collected, full replays passed.
Guarded peak4.0872GB/minfree6.22592GB; charge audit stage observations are separate.
Report outputs/spatial_parts_figure_v2/verification.png was visually inspected;
PDF uses the same figure, not separately inspected. Earlier presentation
variants and every physical source/evidence set remain preserved.
No unit tests were added/run. Full goal and all12 acceptance beads remain open.
The requested earlier 3D viewer and user's original poster remain unchanged.

NEXT read SPECIES_INTERFACE_DIAGNOSTIC_DESIGN.md, then complete the independent
weak species-entropy interface budget, including
one-sided normal-current jumps, potential projection/adjoint and the actual
fluid mass metric. Measure potential/current representation of the fast direction
before any field-refined control; retain explicit B/current transfer errors.
Then pursue the original signed/time/arithmetic, local energy, vacuum, physical
transport/ionization/validity, prescribed-coil/end and converged formation gates.
No solver session needs restarting. The prior status sections below are historical.

## Latest continuation: completed replays and updated field display (2026-10-12)

All physical trajectory, tangent-probe and tangent-audit workers are terminal;
actual exits have been collected. Do not poll or restart session40737 or any
probe worker. The 19-file full trajectory rejection replay, 26-file tangent
probe replay and one-file physical audit replay all match byte identically.
Read TRAJECTORY_TANGENT_VERIFICATION.md. Saved-error spatial projection is
about950741/s, with physical density projection950294/s; these are directional
diagnostics, not certified eigenvalues or an identified physical instability.

At the user's request, the latest evolved fields now have a separate 3D
presentation: outputs/trajectory_display/index.html and current_state_3d.png.
The PNG was inspected. Four real saved times end at28.3333us, accepted step17;
step18 rejected. The display is a filled periodic verification cylinder,
70mm diameter and240mm long, with5mT bias, not a formed FRC or lab design.
Final display-sampled peak|J| is835178A/m² and minimumBz is+4.92882mT.
All four raw display files, metrics and JSON reproduce exactly in a fresh
process. Read TRAJECTORY_DISPLAY.md for independent sampling, interpolation
errors and source records. No physics source or accepted state was changed.

The only live session is the localhost presentation server74537 at127.0.0.1:8786.
It runs no solver. Browser inspection is pending because the Mac was locked;
static PNG visual inspection is complete. The viewer is also standalone offline.
All173 preceding frozen sources and the user's earlier poster remain unchanged.
No unit tests were added or run. Full goal/all12 acceptance beads remain open.

NEXT scientific action: decompose actual spatial-RHS amplification along the
saved early error with unchanged full-sum and wall-reaction checks, then compare
independent continuum and interface terms. Opposite-sign, same-grid arithmetic,
time refinement, full physical/local conservation, heated exact vacuum,
transport/ionization, real coil/end boundaries and converged formation remain.
Do not launch the obsolete seven-job plan unchanged: its first replay is done.
The prior live-worker sections below are historical and superseded here.

## Active trajectory rejection replay (2026-10-12, after c912fdc)

Latest completed work: independent density localization and its fresh full
replay reproduce all216 original density/error fingerprints and all3 derived
raw files bitwise. From13.333us onward more than99.7% of the refined squared
density error is in low axial cell modes abs(k)<=4. Final error spans radius;
first-interface jump is2.27545% in background-normalized surface L2, sampled
maximum4.69502%. Localization does not identify a cause or prove local flux
failure. Initial tiny-error partitions are quadrature-sensitive. Read the new
TRAJECTORY_REFINEMENT_VERIFICATION section and inspected scientific report in
outputs/trajectory_rejection_report. No extra physical trajectory was launched.

After the complete failed-trajectory replay, the proposed next diagnostic is
the full unchanged spatial RHS along the saved early error direction, with
amplitude and independent-mode controls. Read TRAJECTORY_TANGENT_PROBE_DESIGN.md;
its three modules are prepared but have not been executed. Opposite-sign, same-grid arithmetic,
time and complete physical/conservation controls remain required.

The refined F8/64,B32/64 negative100us/60-step attempt rejected step18 after
17 accepted steps (28.333333us). Newton GMRES failed with info=5. The original
stage1e-12/relative1e-10 bounds remain. Session36840/PID89880 is terminal,
exit0 collected through the driver's recorded numerical-failure path. Its
4.439148GiB peak/390.625MiB minimum free disk pass the8GiB/256MiB guards.
This is a failed full-duration gate, not a resource interruption or success.

Independent density-only observations through all17 accepted steps, quadrature
16/20, fully replay every time-series byte and physical fingerprint. At the
last endpoint, the refined density departs16.0233% at a stored node; its L2
error is790.258 times the growing independent linear density signal. Older
medium/coarse controls have relative errors5.700e-5/3.898e-4 at the same time.
Integral quadrature agreement is close; sampled maxima are less resolved.
Single-sign density evidence does not establish causality or full physical
accuracy. Read TRAJECTORY_REFINEMENT_VERIFICATION.md and the companion for
examples.trajectory_density_prefix. Stored terminal initial/final state and
magnetic arrays match the initial/last accepted arrays in dtype, shape and bytes.

LIVE: parent unified session40737/PID93211, replay child PID93212. Parent record:
outputs/trajectory_rejection_campaign/summary.json. It runs only the full
negative rejection replay, with original initial data, equations and strict
checks. Latest snapshot:12 accepted replay steps, not completed verification.
The one-job plan is outputs/trajectory_refinement_diagnostics/rejection_replay_plan.json.
The old seven-job plan is preserved but was NOT launched. The controller and
new density observer are now used/frozen; do not edit their Python sources.

NEXT poll the actual parent40737, without restarting or launching another heavy
job. Collect its terminal exit and verify all initial/final/state/stage files,
non-runtime solver/failure metrics and exact failed step before claiming a
reproduced rejection. Reproduced failure still fails100us. Inspect the result
before selecting opposite-sign, same-grid arithmetic and120-step controls;
retain full requested durations. Do not launch the obsolete seven-job plan as
is, because its first replay is already active. Original preserved medium
controls use older induction arithmetic, so a same-grid split control is needed.
Full signed field/current/E, pressure, conservation/history/energy audits and
mesh/time controls remain open. Resource-budget each later job and preserve
all failures; full remaining campaign storage is not presumed available.

No unit tests. The full goal/all12 acceptance beads remain open. Completed
source and evidence form a local progress checkpoint; evolving live replay
outputs remain outside that completed evidence set.
All168 old sources,251 prior evidence files,504 prior source archives and the
user's poster were rechecked unchanged. Active173-file source:
8dc9d2e9557b53a27351c973e71762c0577837e4b0024bab1c718f7b209066d8.
See localization_preservation.json and active_source_v2 for the new closure.
See primary_terminal_preservation.json, primary_exit.json and live_progress.json
under outputs/trajectory_refinement_diagnostics. Earlier all-workers-terminal
statements below are historical. This is progress, not a blocked or complete goal.

---

# Continue from complete refined amplitude controls

All workers are terminal and collected. No session needs resuming. No unit
tests were added or run; retain the user's end-to-end-only preference.
The full FRC goal and all12 acceptance beads remain open. This continuation
from cd66ec5 made verified progress.

Read the newest GOAL_STATUS, REFINED_CONTROL_DESIGN/VERIFICATION and
outputs/refined_control_diagnostics/checkpoint.json. Scientific PNG/PDF:
outputs/refined_control_report; PNG visually inspected. The old live snapshots
under refined_control_diagnostics are historical; physical_completion.json and
the final checkpoint supersede them. Do not restart those finished jobs.

## Verified results and limits

- Four accepted amplitudes: epsilon1e-5/1e-3 on B256/320 and B320/256, fixed
  F64/256 and independent radial16/N128 reference, original tolerance1e-10.
  Every full signed workflow and fresh separate-process replay passes.
- Both grids at epsilon1e-4/tolerance1e-11 strictly reject on minus after all
  five outer attempts; final dual-Faraday6.14598e-11/8.29280e-11. Complete
  rejection replays match all input arrays and metadata. Saved rejected inputs
  equal the baseline in shape/dtype/bytes, not a post-call immutability proof.
- Four component, two spectrum and two magnetic/current-curl audits replay.
  Axial-grid Jr error1.76199% becomes4.17506% for epsilon1e-5; strong curl
 1.73644% becomes3.19274%. Radial-grid Jr5.90438% becomes6.70354%. Actual
  amplitude differences3.83376%/3.19588% are hidden by small total-J errors.
  Larger amplitudes slightly worsen radial/axial current and total-J errors;
  azimuthal current barely changes. No amplitude-convergence claim.
- Above cell mode|k|64, Jr error grows9.9952x/9.7616x as epsilon shrinks10x.
  That band holds93.7982%/87.2576% of squared amplitude difference. This fits
  numerical noise amplified by differencing but establishes no cause. Lower-
  band errors persist. Cell spectra fold aliases; strong curl omits face jumps.
- All240 spectrum/component and48 current/component-curl integrals match.
  Volume component quadrature differences<=3.440e-15; band differences reach
 4.968e-5. Parseval/curl decomposition defects<=3.150e-16. Instantaneous
  work defect<=3.983e-8/exchanged powers; local and time budgets remain open.
- Maximum physical peak7.574GiB, minimum physical disk605.469MiB satisfy all
 8GiB/256MiB guards including final extrema. Audit resource checks occur at
  stage boundaries only. All shared completed replay paths remain immutable.
- Fourteen replay groups /14 raw comparisons;58 new archives,446 preceding
  archives,164 preceding sources,228 prior evidence files and22 active older
  source closures verify. Four new modules have companions. Final168-file SHA:
 5deb5dc6c4a3229d60d13d33b2a182716aadd61fcd1fcfa1fa5264ef6dff4966.

## Concrete next action: complete trajectory refinement

1. Read JOINT_E2E_DESIGN.md, PHYSICAL_INTERFACE_VERIFICATION.md and companions
   for matched_mode_fixture, matched_mode_evolution, split_mode_evolution and
   split_mode_audit. Used Python sources are frozen; add a new guarded opt-in
   orchestration module if necessary, with explicit dynamic source seeds.
2. Use the unchanged split_mode_evolution full nonlinear protocol on
   elements8,points64,field_nr32,field_nz64, retaining the existing radial1/N64
   reference and original physical settings. Run both signed epsilon1e-4
   perturbations for100us with60 and120 Gauss steps, original stage/endpoint/
   induction bounds and complete saved stage/state records. Preserve failures;
   a partial run cannot pass the full-duration gate. Run one large job at once.
3. Compare to preserved F4/64,B16/64 trajectories
   outputs/matched_mode_guarded_ne4_z64_t60_minus and _plus, keeping physics,
   duration, reference and amplitude fixed. A solver-arithmetic difference must
   be explicitly audited, not called a bitwise cross-implementation match.
   Fresh full replays of each new trajectory, independent q16/20 field/pressure/
   current/electric/energy/particle/entropy audits and time differences are
   required. The older4.31611% final-current error and fast near-grid branch
   remain unresolved; refining is not presumed to repair them.
4. Source-derived retained dense modal storage0.0579GiB and60-step secant
   history0.0397GiB are lower estimates, not peak limits. Keep actual8GiB RSS
   and256MiB disk guards and final-extrema checks; budget all saved states,
   replays and Git evidence. The finest instantaneous grid is not a feasible
   trajectory simply because its signed rates fit memory.
5. Keep refined-grid noise/strict1e-11 rejection and joint B320/320 resource/
   accuracy controls open. Any compensated arithmetic or expanded electric
   representation requires unchanged real equations/adjoints/bounds and full
   physical controls. No filters, floors, artificial mass or repairs. Continue
   local budgets, heated moving exact vacuum, transport/ionization/validity,
   prescribed real coils/ends, converged formation and experimental UI.

Preserve every archive and the user-modified outputs/frc_showcase/poster.png,
SHA a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca.
No remote or publication; local issue export only. Provisional tube/gas/fill/
bias inputs are not validated laboratory recommendations. No3-D stability or
particle-confinement claim. Older handoffs below are historical.

---

# Continue from independently refined field directions

All workers are terminal and collected. No process or session needs resuming.
This continuation from a67d240 made verified progress. No unit tests were added
or run; continue with complete end-to-end physical workflows, independent
comparisons, convergence/conservation measurements and fresh full replays.
The full original FRC sandbox goal and all 12 acceptance beads remain open.

Read newest GOAL_STATUS, SCHEDULED_WORKSPACE_DESIGN/VERIFICATION and
outputs/scheduled_workspace_diagnostics/checkpoint.json. Scientific PNG/PDF:
outputs/directional_workspace_report; PNG visually inspected.

## Verified evidence

- Scheduling alone preserves coarse/fine physical results bitwise but peaks
  6.481/6.476 GiB. A read-only full replay identifies a second transport peak
  during full stress construction. ComponentDriftStress supplies original
  ((m_reduced/e**2)*J_i)*J_j/n entries to the unchanged stress adjoint, retaining
  arithmetic association and independent off-diagonal entries. Complete
  original/new coarse/fine comparisons and all candidate replays are bitwise.
  Fine primary/replay peak 5.342/5.516 GiB; defaults/time integrators unchanged.
- Fixed F64/256, independent radial16/N128 reference, epsilon1e-4, bound1e-10.
  New B256/320 axial and B320/256 radial signed workflows and full fresh replays
  pass every original check and independent q16/20 physical/face audit.
- Axial refinement Jr errors: 5.96939% (256), 3.05317% (288), 1.76199% (320).
  Independent strong-curl Jr at320 is1.73644%; weak/strong difference0.131541%
  slightly grows. Br error falls to0.0266691%. Modest-grid log slopes5.692/5.218
  do not prove asymptotic convergence or resolve amplitude/cancellation error.
- Radial320 instead gives Jz0.734849% versus1.18700%, strong-curl Jz0.599842%
  versus1.00320%, and weak/strong Jz0.210047% versus0.345296%. Jr stays5.90438%.
  Radial-velocity error improves0.183143% to0.0983390%. Total-J0.0552027%
  still hides weaker components. Strong curl excludes face distributions.
- New maximum residuals6.82535e-11/8.96429e-11 stay below1e-10. Largest
  instantaneous power defect4.343e-9/exchanged powers; local/time budgets open.
  New q16/20 component differences<=2.64e-16; physical rate norms<=1.39e-15.
- Axial primary/replay peak6.901/7.162 GiB; radial7.161/7.037 GiB. An earlier
  axial replay hit the disk guard (262144000 free bytes), exit130 collected.
  Partial physical/source/control records are preserved; no solver rejection.
  Fresh standalone replay succeeds. Final successful guards all obey8GiB and
  256MiB limits; do not run a large solver alongside a diagnostic audit.
- Twelve replay groups: six complete signed simulations, four implementation
  comparisons, two independent curl audits. All8 raw comparisons,49 new archives,
  397 preceding archives and22 older active source closures verify. Eight new
  modules have companions. Final164-file source SHA:
  a7022b5b4d1a8f5b598c431e71f6865d51c178ca61eed43a33ef230cfaeeabf5.

## Concrete next action: refined controls, then joint accuracy

1. Read component_stress_transport, guarded_component_stress, current_rate_curl_audit
   companions and previous amplitude/curl/residual verification. All used Python
   files are frozen. Existing guarded_component_stress accepts nr/nz but fixes
   epsilon1e-4 and tolerance1e-10. Add a NEW explicitly sourced opt-in protocol
   for epsilon1e-3/1e-5 and tighter residual controls on B256/320 and B320/256.
   Keep both signed full rates, exact actual inputs, every original residual,
   q16/20 physical/face norms, component differences and full fresh replays.
   The preceding 1e-11 rejection and raw-curl arithmetic evidence remain valid;
   do not loosen a bound or assume refinement repairs that cancellation.
2. Preserve component and actual magnetic-curl norms. Separate amplitude
   differences from continuum error; the small total-J norm can conceal Jr/Jz.
   Keep point/axis/wall and initial face diagnostics distinct from evolved
   local charge/enthalpy/Poynting balances. No filters, atmosphere or repairs.
3. Budget combined B320/320 and further independent refinement only after
   measured resource headroom. Radial field count must subdivide64 fluid cells.
   The unchanged8GiB/256MiB guard must include final observations. A combined
   run's feasibility is unproved. If storage is limiting, profile the component
   variant in a required full replay before choosing another workspace change;
   original UniformWallMetric modal factors grow roughly as nr^2*nz.
4. Continue full signed trajectory mesh/time/residual convergence and local
   budgets. The radial1 F8/64,B32/64 full100us control remains pending. Rate-only
   fidelity does not establish time stability or permit formation adoption.
   Then resume heated moving exact vacuum, transport/ionization/model validity,
   real prescribed coil and physical-end boundaries, converged formation and UI.

## Data and execution

New physical source:000fd344b524568a655485308fb72d754aa8c144d26cf3a8e04c4d947ffb53ef.
Fine baseline: outputs/component_stress_fine; axial: component_axial320;
radial: component_radial320. Successful axial replay is suffixed _replay_serial,
not the preserved interrupted _replay. Each curl audit has _curl, _curl_replay
and _curl_verification directories. Diagnostics include complete source/replay
inventories, directional metrics, interruption and verified storage manifests.

All byte-identical shared NPZ paths are immutable; never overwrite any in place.
A cross-study inventory rehashed1776 completed paths before/after sharing852
additional duplicate inodes. No scientific content/path or user file was lost.
Do not delete partial attempts or retry into an existing destination. All old
source/evidence and the user's modified outputs/frc_showcase/poster.png survive.
No remote or publication. Provisional lab dimensions, gas, neutral fill and bias
remain exploration inputs; no3-D stability or particle-confinement claim.
Older handoffs follow as historical records, not active worker instructions.

---

# Continue from fully replayed axial refinement and workspace verification

All workers are terminal and collected. No handle needs resuming. No unit tests
were added or run; continue with full end-to-end physical workflows, independent
comparisons, convergence, conservation and complete replays. The full original
FRC goal and all12 acceptance beads remain open. This turn and5a60b36 both
made verified progress.

Read newest GOAL_STATUS, WORKSPACE_MEMORY_DESIGN/VERIFICATION, and
outputs/workspace_memory_diagnostics/checkpoint.json. Scientific PNG/PDF:
outputs/workspace_refinement_report_layout; PNG visually inspected.

## Verified evidence

- Releasing expired current-workspace temporaries and detaching completed
  density models from geometry cache templates preserves every original coarse/
  fine input, rate, face array, solver metric and independent q16/20 norm bitwise.
  Four original/new whole-workflow comparisons and fresh candidate replays pass.
  Original functions, default geometry and time integrators remain unchanged.
- Finest baseline remains F64/256,B256/256 on independent radial16,N128,
  epsilon1e-4,tolerance1e-10. The detached-cache baseline and full replay peak
  6.446/7.056GiB; lifetime-only6.873/7.430GiB shows resource variability.
- New F64/256,B256/288 full signed case and exact complete replay pass all
  original checks, worst6.63361e-11 under1e-10. Peaks7.268/7.408GiB obey the
  unchanged8GiB/256MiB guard; no guard fires. Do not infer doubled-grid feasibility.
- Independent current-component and actual magnetic-rate curl audits reproduce
  original fingerprints and replay fully. Jr error5.96939%→3.05317%; strong-curl
  Jr5.94247%→3.03377%. Weak/strong difference0.110758%→0.120095%. Thus the
  actual magnetic response improves, without a repaired weak Ampere identity.
  Jz remains1.11402%; full-J error0.0817008% hides weak components. Br error
  halves to0.0491827%; Btheta improves to0.363014%. Two grids are not an order.
- Refined instantaneous power defect<=3.111e-9 of exchanged powers. Largest
  q16/20 current-norm change4.76e-16 absolute; physical rate norms1.30e-15.
  Strong curl excludes distributional jumps; local/time budgets remain open.
- All11 replay groups and7 raw-file comparisons pass;9 new modules have
  companions;42 new source archives,355 earlier archives and22 older closures
  verify. Final156-file source SHA7e89bbaac2a377a3abe1d2a6fab0b0641381dbd777ae04a903a5c86e8c32e918.

## Concrete next action: memory and directional accuracy

1. Inspect released_current_transport.py/.md, physical_derivatives.py/.md,
   detached_metric_wall_plasma.py/.md and WORKSPACE_MEMORY_DESIGN. The released
   pair still constructs its full velocity-gradient tensor before constructing
   and evaluating the full stress tensor. These independent evaluations retain
   both9-component tensors together. In a NEW opt-in subclass, move gradient
   construction until after stress_rate completes and stress is released. Keep
   each arithmetic expression, state, density, returned term, face contribution,
   physical mass and original check. This changes call scheduling, so do not
   reuse the earlier del-only AST-identity claim. Verify full coarse/fine old/new
   physical workflows and exact replays, and measure actual peak storage.
2. If rescheduling is insufficient, evaluate stress/gradient contractions by
   component or stream exact quadrature batches in another frozen-source variant.
   Preserve all cylindrical connections and force/work adjoints. Any changed
   floating-point reductions require actual physical component/work comparisons,
   not an assumption of equivalence or a relaxed bound. No floors or repairs.
3. Budget axial320+ and independent radial refinement. Radial field cells must
   subdivide the fixed64 fluid elements, so320 is the next admissible radial
   increment; it also increases the dense uniform metric more than axial
   refinement. Retain the improved self-terminating guard pattern and one large
   solver process at a time. The completed256→288 axial step does not replace
   doubled grids. Old B256/512 resource-interrupted attempt stays untouched.
4. Require full signed cases, independent q16/20 physical norms, each current
   component and actual Bt strong curl, face/local budgets and fresh full replay.
   Add amplitude and residual controls on new grids before declaring convergence;
   epsilon1e-5 sensitivity/raw-curl cancellation is not solved by these memory
   changes. Further fixed-field fluid refinement may be needed; current results
   do not separate every source of remaining Jz/radial-fluid error.
5. Continue full signed trajectory mesh/time/residual convergence and local
   budgets, then heated moving exact vacuum, transport/ionization/validity,
   actual prescribed coils and physical ends, converged formation and the UI.
   The radial1F8/64,B32/64 full100us control is still pending; near-grid growth
   in old trajectories is unresolved. Do not adopt new rate-only evidence as
   time stability, confinement or a laboratory experiment recommendation.

## Practical continuation

All current new files are frozen after calculation. Extend them rather than
editing a used Python source. Complete physical source identity for the detached
workflow is f2402996f3cdf2aca408a538ffaac8e850fea0ba964d7fb6ac8d1f385b1c43ae.
Baseline outputs: detached_current_fine; refined: detached_axial288; independent
current audits: detached_baseline_curl and detached_axial288_curl, all under
outputs/. Existing guarded_detached_refinement fixes F64/256 and accepts field
--nr/--nz; --original verifies a fresh same-source replay, not an algorithm change.
Use streamed_audit_comparison for whole-workflow bitwise old/new comparisons.

Two new manifests share496702418 nominal duplicate bytes from completed,
verified, byte-identical archived NPZs. Preserve all paths/hashes and never
modify either shared path in place. Check actual disk headroom before new runs;
no system-wide storage cause is asserted. No unrelated file was deleted.
Preserve outputs/frc_showcase/poster.png (user hash a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca).
Use bd --no-daemon and local sync --flush-only. No remote/publication. Provisional
lab envelope remains~1m quartz tube,10–20cm diameter,selectable gas,0.1–100mTorr
neutral fill and adjustable Helmholtz-like bias. Hardware stays prescribed
waveforms; none of these fixtures establishes3-D stability or confinement.

---

# Continue from verified tensor mass inverse and measured resource limits

All workers are terminal and collected; no tool handle needs resuming. The full
original FRC goal and all 12 acceptance beads remain open. No unit tests added
or run: continue with complete end-to-end physical workflows, independent
comparisons, convergence, conservation and full replays.

Read the newest GOAL_STATUS, FIELD_COMPONENT_REFINEMENT_DESIGN/VERIFICATION,
and outputs/field_component_diagnostics/checkpoint.json. The preceding09b9363
checkpoint and this continuation both made verified progress.

## Current evidence

- Full original/new fine baseline and component comparisons pass with fresh
  complete replays. Tensor inverse changes only factorization evaluation;
  original equations, masses, quadrature, boundaries and checks remain.
  Default geometry/time integrators remain unchanged. Coarse full workflows
  capture identical actual assembled operators; fine physical inputs agree.
- Radial current error persists at 5.96939%, despite tiny vector error.
  Fine original/tensor Jr difference is 1.04097e-5 of its reference norm.
  Worst strict residual 6.4967e-11 under1e-10; power defect <=8.75e-9 of
  exchanged powers. This is an instantaneous verification fixture, not formation.
- Existing fluid32/64 same-field comparison now resolves each current component:
  Jr errors11.6635%/5.96939%, Jz36.1204%/1.18700%. Full audit replays bitwise.
- Original B256/512 initialization ended with resource interruption, no signed
  input/rate. Worker36206 exit143 was collected. Never resume/restart from its
  partial summary alone or relabel that preserved run as a numerical failure.
- Tensor full fine baseline/replay peaks7.420/7.150GiB; observation replay7.432.
  CurrentTransport.pair first-sign cumulative peak4.829→6.558GiB; second signed
  electromagnetic pair7.061→7.432GiB. Nested high-water intervals are not
  exclusive/live storage. See tensor_resource_observation/calls.json and full
  verified replay. All successful tensor guards remained untriggered.
- Storage estimates3.094GiB forB256/512 and5.098GiB forB512/256 exclude large
  factor/workspace peaks. Do not infer doubled-grid feasibility from them.
- All37 new archives,318 preceding archives and22 older closures are verified;
  147-file source SHA8ae692975fec4c2fac1345f723b13bb428c1bb1affa2c72557f527d476c8b332.
  All historical calculation Python files and the user's poster are unchanged.

## Next actions: frc-7mm.6 and frc-7mm.7

1. Reduce actual temporary quadrature storage in a new opt-in module. Inspect
   current_transport.py and regular_current.py companions and the source before
   modifying anything. CurrentTransport.pair retains full velocity-gradient and
   stress tensors beyond their last use; advection also forms full gradient
   arrays. The second signed electromagnetic evaluation sets the observed
   global peak, so transport alone is not established as sufficient. Profile
   before/after using own-worker resources and unchanged physical checks.
   Do not edit frozen calculation modules; add a separate variant.
2. Require complete original/new signed coarse and fine workflows, identical
   physical inputs/operators where applicable, independent physical components,
   work, q16/20 and fresh full replays. Arithmetic changes need measured actual
   differences, not assumed bitwise equivalence. Keep all original failures,
   physical masses/density, basis modes, bounds and numerical budgets.
3. Retry fixed-reference F64/256 separate axial B256/512 then radial B512/256
   only with justified memory/disk headroom. Use the new self-terminating guard
   pattern; the old SIGINT-only guard failed to promptly stop native work.
   Require complete both-sign rates, component/strong-curl comparisons and full
   fresh replays. Never substitute total-current norm for Jr accuracy.
4. Continue full signed trajectory mesh/time/residual convergence and local
   budgets, then heated moving exact vacuum, transport/ionization, real prescribed
   coils/ends, converged formation and experimental UI. The radial1F8/64,B32/64
   full100us control remains pending; old near-grid growth is unresolved.

Shared-storage manifests cover only byte-identical completed archived replay
NPZs. Both paths now share immutable bytes: never modify either in place. No
research evidence was deleted. Check actual free disk before new allocations;
host-wide changes are not attributed to a particular process without evidence.
Preserve outputs/frc_showcase/poster.png (user hash a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca).
Use bd --no-daemon and local sync --flush-only. No remote/publication. Provisional
lab envelope remains~1m tube,10–20cm diameter, selectable gas,0.1–100mTorr neutral
fill and adjustable Helmholtz-like bias; hardware remains prescribed waveforms.
These fixtures establish neither 3-D stability nor confinement.

---

# Continue from separated curl arithmetic and current-component error

Read the newest GOAL_STATUS and CURL_EVALUATION_DESIGN/VERIFICATION.
**All workers are terminal and collected. No job needs resuming.**
Use complete workflows, independent physical comparisons, convergence,
conservation and full replays. Do not create or run unit tests. The complete
original FRC sandbox goal and all 12 acceptance beads remain open.

## Verified current evidence

- Original strict-rejection vectors reproduce bitwise. Simple projected-curl
  factoring is not a repair. Exact raw products/sums plus addition-low accounting
  leave 2.97e-15/RHS distributivity error, but actual dual Faraday remains
  5.62e-11 against 1e-11 on that saved iterate. No solver change is adopted.
- Full physical component spectra at epsilon 1e-3/1e-4/1e-5 reproduce original
  fields and replay completely. Jr errors 5.9751%/5.9694%/6.6721%; high-band
  norm errors 0.02750%/0.27966%/2.78699% approximately scale as 1/epsilon.
  The dominant |k|20 spatial-error hump persists. Above |k|64 contains 86.672%
  of the 1e-4 versus 1e-5 amplitude-difference energy, but only 0.2195% of
  the 1e-4 squared error. No accepted state was filtered. These are cell modes,
  not unaliased continuous harmonics. Full original integrals/Parseval verified.
- Differentiating the actual centered magnetic rate independently gives 5.9425%
  radial-current error versus 5.9694% in Jt, with 0.1108% weak/strong difference.
  This is an inside-cell volume comparison excluding distributional face jumps.
  Signed cross terms close the error identity to 2.16e-16. Both amplitudes replay.
- New complete signed epsilon 1e-3 F64/256,B256/256 rates pass every original
  check, max 1.91046e-11 under 1e-10. Full q16/20 physical audits and exact fresh
  replay pass. Increasing amplitude slightly worsens velocity/toroidal-field
  errors: numerical-noise reduction is not an accuracy or formation gate.

## Next actions: frc-7mm.6 and frc-7mm.7

1. Budget independent field-axial and field-radial refinement at the fixed
   radial16,N128 reference and fixed fluid representation. Measure each current
   component, actual Bt derivatives and local budgets. Do not adopt the tiny
   total-current norm as evidence for Jr accuracy. Use the epsilon 1e-4 baseline
   alongside amplitude controls; 1e-3 visibly changes nonlinear responses.
   Finest-mesh memory and sparse-factor storage need explicit accounting.
2. Investigate accurate raw-curl and electric-sum representation. Applying C
   then P alone did not solve the original rejection. Keep original physical
   masses/density, real equations, adjoints/work, kernel terms, all residuals
   and bounds. A proposed arithmetic implementation needs complete exact-input,
   physical-power/field comparisons and signed full replays before adoption.
3. Continue complete signed radial1 F8/64,B32/64 100us trajectory controls
   against preserved F4/64,B16/64 when accuracy gates justify adoption, including
   time/mesh/residual and independent local conservation. Near-grid growth is
   unresolved; do not treat new instantaneous rate audits as time convergence.
4. Continue heated moving exact vacuum, transport/ionization and validity,
   physical prescribed coils/ends, converged formation and experimental UI.
   The user's ~1 m tube, 10-20cm diameter, variable gas, 0.1-100mTorr neutral fill
   and adjustable bias remain provisional. Hardware stays abstract waveforms.

## Evidence and preservation

- outputs/curl_evaluation_probe and its _replay/_verification: all four original
  iterates, alternative full actions and final exact raw reference.
- outputs/component_axial_spectrum and its _larger variant, both with complete
  _replay/_verification: full physical component spectra and amplitude effects.
- outputs/current_rate_curl_comparison and _replay/_verification: independent
  magnetic/current-rate decomposition, q16/20 and actual field fingerprints.
- outputs/radial16_amplitude_larger and _replay/_verification: new full signed
  physical rate workflow. Every coefficient/face array and metric replays.
- outputs/curl_spectrum_report_layout/verification.png: visually inspected.
  Original report layout also preserved. No PDF generated in this continuation.
- outputs/curl_evaluation_diagnostics/checkpoint.json: complete source, input,
  worker and replay inventory. No disk guard triggered. Large geometry builds
  and replays used separate processes.

Five verified groups comprise four diagnostics and one full signed rate
calculation. Six new Python modules are frozen, each with its companion.
18 new archives, 300 earlier archives and 22 earlier closures verified.
Final 140-file source SHA b0ee5f17e5baf0fa4fc564d3fe40437173b05479d16210e2ef040f0a1d0ce922.
No old calculation source changed. Preserve user poster SHA
 a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca
outside research commits. All work remains local; no remote/publication.
Older handoffs below are historical.

---

# Continue from component accuracy and split-residual path diagnosis

Read the newest GOAL_STATUS and SPLIT_RESIDUAL_PATH_DESIGN/VERIFICATION.
**All workers are terminal and collected; no job needs resuming.**
No unit tests: use complete workflows, independent physical references,
convergence/conservation measurements and full replays. Keep all 12 acceptance
beads and the overarching FRC goal open.

## Current evidence

- The original finest tighter-bound failure is reproduced with unchanged input,
  equations, residuals and budgets. New observations show final opposing-vector
  gap 6.87234e-11 of RHS. Electric curl distributivity leaves only 0.4245% of that
  gap; H/T and current/inertia paths do not dominate. All-term remainder is
  5.41e-17 of RHS. The complete fresh standalone observation replay is bitwise.
- Exact Decimal sums verify 786432 actual electric coefficients, with zero
  mismatches and a complete postprocessing replay. Propagating addition loss
  alone leaves 66.21% of the gap unexplained. Do not conflate addition rounding
  with the full curl-action discrepancy or infer an accepted solver repair.
- Complete independent physical component audit/replay on F64/256, B256/256,
  fixed radial16, N128 reproduces every original field/reference fingerprint.
  Jr-rate errors are 0.05969389 / 0.06672077 at epsilon 1e-4 / 1e-5. Jz errors
  are 0.01187001 / 0.01193623; Jtheta error is 0.0001904658. The total J error
  0.000878993 therefore does not certify these smaller physical components.
- Radial-current amplitude difference is 0.03014823 of its own continuum norm.
  88.674% of its squared difference is interior and 11.115% in the outer 5% R.
  Radial-current error itself is about 6% in each band. Axial-current and radial-
  velocity squared errors are more wall-localized (58.08% and 72.05% there).
  All nine fluid/magnetic and all E/J components are measured; q16/20 norm and
  masked-band differences are recorded. Band integrals mask original nodes,
  rather than subdividing cells at the two fixed physical band edges.

## Next actions: frc-7mm.6 and frc-7mm.7

1. Inspect original `StrongWallGeometry.faraday = projection @ original.faraday`.
   On the exact captured original iterates, compare assembled action with raw
   curl followed by the same projection, and isolate product/sum rounding.
   Keep actual electric terms and the independently verified addition low part.
   Do not replace a rejected full rate with a diagnostic vector. Any proposed
   arithmetic change must independently verify the original full equations,
   physical fields and work, preserve strict checks/budgets, and fully replay
   complete signed workflows. The observations are not a repair.
2. Determine whether component errors arise from representation, axial/radial
   differentiation, nonlinear amplitude effects or arithmetic. Use physical
   component norms and fixed independent radial16, N128 reference. Axial-mode
   diagnostics or separate directional refinement can guide the next grid;
   do not infer a cause merely from the dominant total-current norm. A causal
   link to the strict residual-path mismatch remains unproven. Any Fourier
   analysis is diagnostic only and must not filter accepted states.
3. Continue complete signed radial1 trajectories at F8/64, B32/64 against the
   preserved F4/64, B16/64 runs when the accuracy gates justify adoption. Require
   100 us, time/mesh/residual controls, local budgets and independent physical
   conservation. The prior near-grid growing branch remains unresolved.
4. Continue heated moving exact vacuum, defensible transport/ionization and
   validity diagnostics, real prescribed coils/ends, converged formation and
   the experiment UI. Keep original lab inputs provisional and distinguish
   reversal, closed flux, equilibrium, 3-D stability and confinement.

## Artifacts and resource constraints

- `outputs/split_residual_probe`: original full exact-input observation.
- `outputs/split_residual_probe_replay_v2` and `_verification`: passing complete
  standalone replay. The older `_replay` directory is partial and preserved.
- `outputs/split_residual_arithmetic_audit` and component_amplitude_localization:
  original audits plus their `_replay` and `_verification` directories.
- `outputs/split_precision_report/verification.png`: visually inspected
  scientific figure, with a PDF alongside. PDF not separately inspected.
- `outputs/split_residual_diagnostics/checkpoint.json`: source/input/worker
  inventory, resource guards and preserved interruption.

Run large original/replay geometry builds in separate processes. The earlier
in-process replay was interrupted when free disk fell below 256 MiB; no global
cause is asserted and no unrelated files were deleted. The standalone replay
completed with the same physical calculation and guard. Source-derived finest
trajectory modal storage still exceeds 15 GiB before field/work arrays.

All six new modules and old calculation sources are frozen; read companions
and extend via separate modules. Final 134-file source union:
150dc2f3e7614eed5da58d314b6ff3e8cf804b641c1f9fec9f91b738e4712745.
Three complete replay groups and 12 new archives verified; 288 preceding
archives and 22 older source closures unchanged. User poster SHA
 a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca.
Preserve it outside research commits. Beads and research remain local; no
remote or publication. Older handoffs below are historical.

---

# Continue from complete joint refinement and physical amplitude sensitivity

Read the newest GOAL_STATUS, JOINT_E2E_DESIGN and JOINT_E2E_VERIFICATION.
User policy: no unit tests; use complete workflows, independent physical
comparisons, convergence/conservation measurements and full replays.
**All workers are terminal and collected. No job needs resuming.**
All 12 full-goal beads remain open, especially frc-7mm.6 and frc-7mm.7.

## Verified evidence

- Fixed radial16, N128 reference SHA
  3b142fe986aca23d340b1739c17ec83b5b94e04574cff2a661fe1b0a8527ef08.
  Full F32/256 and F64/256 signed cases on B256/256 and full reruns pass.
  Finest physical rate errors: vr 0.00183142646, Br 0.00097962351,
  Btheta 0.00547509477, Jt 0.00087899301, E 0.00000456807405.
- Finest epsilon 1e-5 success and complete replay give Jt 0.00088543967.
  The centered current coefficient difference is 0.00947096, while actual
  physical vector difference is 0.000107450107. Radial/axial current-component
  physical differences are 0.0301482272 / 0.00126286432 relative to their own
  references. E difference is 7.148893e-8. The independent reconstruction
  reproduces original physical fingerprints and fully replays at q16/20.
- The tighter 1e-11 case rejects the identical minus input and reproduces its
  rejection bitwise. No failed rate exists. Saved raw input:
  outputs/radial16_finest_controls/case_01.npz, SHA
  ce34210139b000f4c2746dfde664ba1dd188d5e2485572a364e184bf12c4d999.
  Dual Faraday residual stays near 5.89e-11; scalar norms imply residual-path
  mismatch at least 2.499e-11 at the recorded final iterates. This does not
  identify the cause or prove impossibility at other iterates.
- Both separate radial1, N64 signed 100 us / 60-step trajectories, every saved
  stage/state and all physical/conservation audits replay exactly. Actual
  energy drift is at most 3.636e-7 of initial physical signal-energy scale.
  These coarse integrations do not establish refined time/mesh convergence
  or remove the previously observed near-grid growing branch.
- 20 complete replay groups (including postprocessing, rejected cases and
  completed-case selection), 137 raw-file comparisons and 71 source archives
  verified. The preceding 42/23/59/61/32 archives and 22 older source closures
  remain unchanged. See outputs/radial16_joint_diagnostics/checkpoint.json.

## Next concrete actions

1. Reconstruct the exact rejected input with its current split-inverse and
   direct-metric setup, enforcing bitwise input identity. Capture the actual
   operator/dual residual vectors and individual mass-inverse actions without
   changing the arithmetic or acceptance. An observation-only copy must
   reproduce the original first-attempt rejection exactly before interpreting
   its decomposition. The old balanced WallRoundoffCore is different arithmetic
   and cannot be used as a drop-in observation of this failure. Preserve every
   original bound, physical mass, actual density and refinement budget.
2. Investigate the radial-current amplitude sensitivity in physical component
   norms, including near-axis/wall structure, and reconstruct fluid/magnetic
   difference fields as needed. The existing current/E audit measures these
   fields directly; close total error norms alone do not certify sensitivity.
3. Continue full signed radial1 trajectories on F8/64, B32/64 against preserved
   F4/64, B16/64 results, keeping the independent reference and physical inputs
   fixed. Require full 100 us runs, physical field/conservation audits, complete
   replays and time controls. Preserve failures and growth. This candidate grid
   is not presumed accurate. Budget memory first: source-derived retained dense
   modal storage is 0.0579 GiB here, but 15.414 GiB on F64/256, B256/256, before
   field matrices, sparse factors, nonlinear/Krylov and audit workspaces.
4. After those gates, continue local current/enthalpy/Poynting budgets, heated
   moving exact vacuum, transport/ionization, actual coil/end conditions,
   converged formation and the experiment UI. Keep the full FRC goal open.

## Preservation and local state

All calculation code used here is now frozen, including the 14 new opt-in
modules listed in source_preservation.json. Read companions before extending
it; use new modules rather than editing archived sources. The final 128-file
source union is archived under radial16_joint_diagnostics/final_source.
Earlier interrupted batches and partial replays stay intact; completed-case
selection explicitly preserves the incomplete parent status. The finest
control disk guard did not trigger. No unit tests were added or run.

Scientific figure: outputs/joint_e2e_report/verification.png, visually inspected;
PDF is alongside it. User poster SHA
 a783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca.
Do not include that preexisting poster change in research commits. Work and
beads are local with no configured remote or publication. The provisional lab
envelope is not a selected experiment; no 3-D stability or confinement claim.
Older entries below are historical.

---

# Continue from directional term accuracy and strict arithmetic recovery

The full FRC experiment-sandbox goal and all12 acceptance beads remain open.
**470 tests passed in233.11s. Every worker is terminal and collected.** Read the
newest GOAL_STATUS.md, DIRECTIONAL_TERM_DESIGN.md and
DIRECTIONAL_TERM_VERIFICATION.md before changing code. No local or remote job
needs continuation or restart.

## Current evidence

- Same physical radial16,N128 reference, growth12307.1561445/s; original reference
  array SHA3b142fe986aca23d340b1739c17ec83b5b94e04574cff2a661fe1b0a8527ef08.
  Axial doubling improves Br .051938 to .00231942 but leaves J .0847343.
- Radial field128->256 at fluid32/128 lowers vr .1088365->.0309084 and
  J .0875688->.0337504. Br .051625 and Btheta .279381 remain large.
  Fluid32->64 alone improves pressure error31.53-fold but leaves magnetic-force
  error. Refine the two approximations independently rather than attributing
  the full error to one term from a single net residual.
- Independent original pressure/Lorentz and ideal/Hall/pressure/inertia terms:
  full captured fluid sum bitwise matches saved rates; physical terms/curls
  checked independently. Hall/pressure toroidal norms ~24000 times net response.
  Inertia includes numerical-current feedback. Raw errors, cross integrals,
  strong/projected curls, fixed physical wall bands and q16/20 norms remain.
- Original radial256 solve failed with operator2.60e-10 and dual2.83e-10 under
  1e-10. Original failure and exact input77e43521... remain preserved. HT
  association was ruled out; mass-inverse distributivity dominates opposing
  residuals. Independent reconstruction retains7.06% vector remainder.
- Opt-in SplitInverseWallPlasma/Core uses E=M^-1f+M^-1(IJt), sharing forcing
  inverse with original RHS. Every original equation, physical mass, actual
  density, boundary, operator, bound and budget remains. Full exact input passes
  at5.27777e-11, bitwise replayed. Polynomial full B/J/E/work and atomic tests pass.
  Full power defect3.99063e-9 and restored full Faraday rate-relative2.722e-7
  remain explicit; solving accepted equations does not certify physical accuracy.
- Three full signed split-inverse pairs replay bitwise: baseline, radial256
  eps1e-4, radial256 eps1e-5. Physical errors nearly agree across amplitudes;
  centered current coefficient difference .00108516 is retained. The old/new
  baseline inputs and full fluid rates match bitwise; J coefficient change
  5.86e-10. The guarded trajectory sourcef255a576 remains unchanged/unadopted.
- Thirteen complete replay groups,42 new source archives verified;23/59/61/32
  preceding archives byte-unchanged and22 earlier source closures unchanged.
  Full checkpoint: outputs/radial16_directional_diagnostics/checkpoint.json.
  Visually inspected report: outputs/radial16_directional_report_v2 PNG/PDF.

## Immediate next action: frc-7mm.7 and frc-7mm.6

1. Read companions for split_inverse_wall_core/induction/plasma,
   split_inverse_radial_rates, radial_parent_terms, radial_term_capture/audit,
   radial_drive_capture/audit and radial_error_localization. All now have frozen
   source closures; do not edit code imported by archived or live calculations.
2. Run a fresh combined-grid comparison F32/256,B256/256 and F64/256,B256/256,
   retaining fixed independent radial16,N128 reference and physical masses.
   Use both signed amplitudes1e-4 and1e-5 where possible, q16/20, strict original
   residuals and explicit failed records; full fresh replay is required. Separate
   initial representation, complete rates, pressure/Lorentz cancellation,
   current and one-sided face errors. Never make a rejected state look accepted.
3. Localize any remaining pressure/current/wall error before adopting the new
   arithmetic into complete growing-response trajectories. Then repeat actual
   signed100us mesh/time/residual convergence and continuous conservation/local
   enthalpy/Poynting/charge face budgets. Old time order remains unestablished.
4. Continue heated moving exact vacuum, transport/ionization and validity,
   physical ends/coils, converged formation and experimental controls/UI.
   No atmosphere, filtering, artificial mass, heat or energy repair.

The roughly1m quartz tube,10-20cm diameter, selectable gas,0.1-100mTorr neutral
fill and adjustable Helmholtz-like bias remain provisional lab variables.
The current smaller argon filled-periodic case is a numerical verification
fixture, not an experiment recommendation. Axisymmetric closed flux does not
establish three-dimensional stability or particle confinement.

No remote is configured; all research and issue state stay local. Preserve
outputs/frc_showcase/poster.png (SHAa783f7cda0f4a249cc4c23ead61debbdf89d356436a2a3dab2189029c6d5deca).

---

# Continue from fixed radial16 accuracy and strict solve recovery

The complete FRC experiment-sandbox goal remains open. Read GOAL_STATUS.md,
RADIAL_MODE_DESIGN.md and RADIAL_MODE_VERIFICATION.md first. **462 tests pass
in 230.35 s.** The active new entry point is examples.balanced_radial_rates.run
for further instantaneous studies; do not silently change guarded trajectories.
Older handoff sections below are historical, including their old next actions.

All jobs from this continuation are terminal and collected. Seven full bitwise
verification groups and23 new source archives are verified. The59 preceding,
61 older matched and32 older wall archives are byte-unchanged, and22 earlier
verified source closures still match. No worker needs resuming. The full
checkpoint is outputs/radial16_diagnostics/checkpoint.json.

## Current scientific evidence

- Fixed independent radial16 reference at N128, checked against N256. Growth
  12307.15614451037/s; full physical/J/E component shape differences <6.9e-9.
  Raw reference outputs/radial16_reference/continuum_01.npz SHA256
  3b142fe986aca23d340b1739c17ec83b5b94e04574cff2a661fe1b0a8527ef08.
- The original N64 archive remains bitwise unchanged but has a .003251 off-grid
  curl/current mismatch from its real-Nyquist convention. Do not use it as a
  resolved high-radial reference. New N128 independent curl checks pass.
- Separate fluid/field/axial meshes show large full-RHS refinement effects.
  F16/128,B64/128 vr/Br/Jt errors3.86912/.129238/5.10160 improve to
  .1088365/.0519380/.0875688 at F32/128,B128/128. **Btheta error .383069**
  remains; this is not an accepted spatial accuracy or formation gate.
  Wall-region vr/Btheta/Jt errors .576775/1.154478/.438075 identify a concrete
  next diagnostic target; these normalize to the local reference norm.
- Initial final n/vr errors~1.05e-5, Bz .000116518, weak J .000266164;
  radial normal-J/tangential-B jumps .000320212/.000323249. Axial jumps
  remain9.93422e-5/.000126339. Strong Ampere .000982571 uses strong J norm,
  excludes distributions and differs from the old extracted-mode denominator.
- The original F32/B128 minus solve strictly rejects at operator1.39069e-10
  versus1e-10. Its raw input SHA256 is
  e7a3cf637ceff0f3201a88f12ddd84b593f28d23f898406915c8d040a3eac755.
  BalancedWallCore proposes the mean of operator/physical dual corrections
  when the operator bound fails. Both original residual paths then pass;
  worst7.72779e-11. All physics, tolerances and refinement budgets unchanged.
  Exact saved-input recovery and baseline bitwise equivalence are recorded.
- Independent analytic cylindrical B/J/E and work tests pass on6/12/24 cells;
  iteration failure is atomic. Full actual power cancellation on the saved
  input is2.696e-9 relative to absolute terms, separately measured.
- Amplitude1e-4/1e-5 gives nearly the same physical error norms; centered
  coefficient Jt differs .000773103 relatively. Keep that limitation explicit.
  Volume quadrature16/20 norms differ <=3.05e-11 relatively. The largest
  axis relative discrepancy belongs to7.84e-11 error and is1.12e-15 absolute.
- Inspected scientific PNG/PDF: outputs/radial16_report/verification.
  Full inventory: outputs/radial16_diagnostics/checkpoint.json.

## Concrete next actions

1. Read companions for archived_contact_mode, radial_contact_fixture,
   balanced_radial_rates, radial_mode_interfaces, balanced_wall_core,
   balanced_wall_induction and balanced_wall_plasma. These are opt-in; all
   earlier source closures and current accepted trajectories remain frozen.
2. Hold the same radial16,N128 reference fixed. Next compare axial refinement
   dict(elements=32,points=256,field_nr=128,field_nz=256,epsilon=1e-4)
   separately from radial-field refinement
   dict(elements=32,points=128,field_nr=256,field_nz=128,epsilon=1e-4), then
   refine fluid elements to64. Account for quadrature20 memory before launching
   concurrent large cases. Every destination must be new and immutable.
3. Localize remaining vr/Btheta/Jt rate error by original pressure/current/work
   contributions and physical wall/interface regions. Retain actual inputs,
   signed rates, full E/J, independent quadrature, amplitude, strong/weak/face
   errors and exact rejection. Do not fit fields or remove species feedback.
4. Once spatial accuracy supports it, repeat complete signed100us response
   mesh/time/residual controls with local budgets and small-signal energy scale.
   The old evolution source is stillf255a576 and has not adopted balanced
   induction. Previous medium final J error .0431611 and unestablished time
   order remain active limitations; rate tests do not supersede them.
5. Continue full heated moving exact vacuum, transport/ionization, real coil/end
   boundaries, prescribed-waveform formation and experiment UI. All12 unfinished
   beads, including frc-goal, remain open until actual acceptance evidence.

Use the repository .venv with PYTHONDONTWRITEBYTECODE=1,
VECLIB_MAXIMUM_THREADS=1, MPLCONFIGDIR=/private/tmp/plasma-bench-mpl,
MPLBACKEND=Agg. New rates use examples.balanced_radial_rates.run; repeat from
scratch and verify with examples.wall_diagnostic_replay.verify. Do not modify
any imported source of a running or frozen calculation. Use bd --no-daemon;
sync issue state locally. No remote is configured; no publication is authorized.

Preserve the user's modified outputs/frc_showcase/poster.png. The proposed lab
remains~1m quartz tube,10–20cm diameter,0.1–100mTorr neutral fill, selectable gas
and adjustable Helmholtz-like bias. No lab design,3-D stability or confinement
has been established. This remains useful ongoing research, not a blocked goal.

---

# Continue from the radial near-grid response failure

The complete FRC experiment-sandbox goal remains open. Read the newest
GOAL_STATUS.md, PHYSICAL_INTERFACE_AUDIT_DESIGN.md and
PHYSICAL_INTERFACE_VERIFICATION.md. **457 tests pass in254.80s.** No production
physics/time-solver code was edited. New source closures are opt-in diagnostics.

## All jobs terminal

Both larger-grid60-step/100us primary runs and complete reruns finished; every
accepted state/stage and non-runtime metric compares bitwise. Their verification
records are outputs/matched_mode_guarded_ne4_z64_t60_plus_verification and the
corresponding minus_verification. Former worker sessions10675/67242 exited0 and
were collected; do not restart or poll them. The coarse120-step pair, both full
reruns, every physical audit, history, interface, extracted-mode diagnostic and
regression are also terminal. No job was killed or accepted with relaxed bounds.
All59 new source archives are verified;61 preceding and32 older archives plus22
earlier verified source closures remain unchanged. See the full checkpoint in
outputs/physical_interface_diagnostics/checkpoint.json.

## New evidence

- Four additional complete100us signed runs, coarse120 steps and medium60steps.
  Medium F4/64,B16/64 intended projection improves to2.40217, but final J error
  is.0431611 and Bz.00496073. This is a failed accuracy gate; no mesh convergence.
- Coarse30/60/120 controls do not establish fourth-order time convergence.
  Representative uniform spectra have omega*dt118/318 at60steps on coarse/
  medium grids. Phase underresolution is a risk, not the sole diagnosed cause.
- Full physical histories replay every saved checkpoint atq12/16, matching the
  separate finalq16 audit. Medium extra current grows~71682/s after cancellation
  of the independent slow reference; direct actual spatial RHS projects~71426/s.
  Probe amplitude1e-4/1e-5 current rate difference1.57e-6; q16/20 agrees.
- Extracted current is97.83% azimuthal, B98.07% axial; axial Fourier1 dominates.
  Radial Bessel15/16 dominate density/Bz on16 field cells.32 modes capture~94%
  of density/Bz. Weak J versus strong volume curl(B)/mu0 differs.599823; this
  excludes distributional face terms. Independent matched continuum action has
  large radial/induction discrepancies; no validated physical eigenmode.
- Exact fluid subdivision to8/16 elements preserves fields<7e-14 and leaves
  current rate nearly unchanged, but changes momentum rate dramatically. It is
  a separate diagnostic, not an established repair. All three extracted-mode
  studies fully replay every array and non-runtime metric.
- Eight actual original physical face/volume budget cases and axis/wall audits
  replay bitwise atq16/20. Energy cell L2 defect.60485 ->.00017354; final
  particle/current/eentropy/ientropy defects1.306e-7/3.791e-7/9.741e-5/4.022e-8.
  Energy absolute trace jumps5.19e-4 ->1.84e-7 W. Max quadrature norm discrepancy
  9.64e-10. These are rate consistency results, not full nonlinear acceptance.
- Initial conservative energy Hessian is independently verified against70digit
  original EOS differentiation and analytic derivatives. Actual total-energy
  drift/initial perturbation scale<=6.42e-7 across five signed pairs; the scale
  is not an independently conserved free energy. All energy audits replay.
- Inspected report: outputs/physical_interface_report/verification.png and PDF.
  New full regression: outputs/physical_interface_diagnostics/pytest_final.log.

## Required next work

1. Resolve the near-grid radial branch using independently matched fluid/field
   radial spaces and fixed-physical-wavenumber refinement, retaining exact
   physical masses, full species/current/work coupling and all wall constraints.
   Concrete first control: the independent radial16,N64 mode is already frozen
   in outputs/stationary_contact_audit/continuum_07.npz (SHA256
   f9c92d8d94be18b10955fdde1f20868e2d3f0fafcc58a7157eca3e318258f6c6),
   with growth12307.156241/s. Add an opt-in reference/fixture that accepts this
   archived mode without changing continuum_contact_mode or matched_mode_fixture.
   Compare actual full rates on field radial16/32/64 (axial64 then128), varying
   fluid elements4/8/16 separately and keeping the same physical radial16 mode.
   Preserve exact initial representation, full E/J, strong and face errors and
   epsilon1e-4/1e-5 controls; refine the independent axial reference separately.
   A comparison that simply follows a different largest eigenvalue as the grid
   changes is insufficient. Include strong/weak/face consistency and uncaptured
   norms; no fitted eigenmode, filter, artificial mass or energy repair.
2. Repeat complete growing signed trajectories with genuine mesh/time/amplitude/
   residual convergence. Keep small-signal energy and actual local interfaces.
   Source for current guarded evolution remainsf255a576; extracted diagnostics
   do not change the accepted trajectory solver or physical model.
3. Continue heated moving exact vacuum (never an atmosphere), transport and
   ionization validity, real coil/end boundaries, prescribed-waveform formation,
   experiment controls/UI. Keep frc-goal and all current acceptance beads open.

Preserve outputs/frc_showcase/poster.png, which was already user-modified. No
remote is configured and no research data was published. The user's~1m quartz
tube,10–20cm diameter,0.1–100mTorr neutral fill, selectable gas and adjustable
Helmholtz-like bias remain exploration inputs. No laboratory formation,
three-dimensional stability or particle-confinement claim follows.

---

# Continue from matched growing-response verification

The full user goal remains active. Read the newest GOAL_STATUS.md,
MATCHED_MODE_SOLVER.md and MATCHED_MODE_VERIFICATION.md. **446 tests pass in
234.86 s.** All new jobs are terminal. All six signed trajectories replay every
array and non-runtime metric bitwise. See outputs/matched_mode_diagnostics/checkpoint.json
for archive and final verification details. No running job needs continuation.
The overall scope has not changed.

## New evidence

- Eight independent field-direction rate cases, all arrays/physical fingerprints
  and both quadratures bitwise replayed. Final Br/Jt errors0.000244246/0.000219551
  at field64/128, fluid16/128; source 7bb77514. The old0.55/0.48 gate is resolved
  as an instantaneous accuracy study, not full fine-grid trajectory validation.
- Modal uniform physical block inverse and physical uniform current-inertia
  preconditioner preserve all actual equations, full Fourier modes and strict
  bounds. Focused independent operator/stage/atomic tests and full regression.
- Six full 100 us signed runs on coarse fluid 4/32 and field 8/32: 30/60 steps at
  epsilon 1e-4 plus30 steps at2e-4. Source f255a576. Three complete independent
  physical audits and two physical time/amplitude comparisons replay bitwise.
- At60 steps n/Br/J/E errors 0.001011/0.001937/0.002748/0.008678. Growth projections
  near 2.402 agree with independent 2.402175837, but initial representation is
  explicit and this alone is not convergence or a confinement result.
- Half-dt changes E by 0.00136156; double amplitude changes all measured responses
  <=2.07e-6. Energy normalized by total background is near roundoff, not a
  small-signal conservation proof. Full physical quadrature 16/20 norms agree.
- GuardedModalGauss bounds unproductive proposal work without changing acceptance.
  Prior slower pilot algorithms and one failed negative trajectory are preserved.
  Eleven pilot accepted stages/endpoints have fresh independently replayed audits.
  A stale-import postprocessor failed before producing physical audit output;
  its frozen on-disk source limitation is explicit in its termination.json.
- Report PNG/PDF visually inspected. All 22 old verified calculation closures and
  all 32 preceding checkpoint archives unchanged. Preserve user's modified
  outputs/frc_showcase/poster.png. No remote or publication.

## Immediate next action — frc-7mm.7 remains open

1. Complete full response mesh/time refinement. Start with
   examples.guarded_mode_evolution.run on spec
   dict(elements=4,points=64,field_nr=16,field_nz=64), both epsilon +/-1e-4,
   duration1e-4 and60 steps. Then independent matched_mode_audit and full
   matched_mode_replay. Every step has an immutable accepted checkpoint.
   Finer target: field32/128, matching fluid axial128, with radial element
   refinement. Fine rate results do not replace actual nonlinear trajectories.
2. Add a 120-step coarse pair to quantify temporal convergence; use the physical
   comparison, not merely a fitted amplitude. Test tighter residuals separately
   and preserve any rejection. The original extra-tight1e-14 prototype rejected;
   never silently relax bounds or claim it passed. Measure initial small-signal
   energy scaling before interpreting background-normalized roundoff drifts.
3. Close actual local current/enthalpy/Poynting interface and jump budgets.
   Continue heated moving exact vacuum (no numerical atmosphere), transport/
   ionization and model validity, physical ends/coils, converged formation and
   experimental controls/UI. Keep all acceptance beads and frc-goal open.

All analysis is fully ionized bounded periodic verification. The user's roughly
1m quartz tube, 10–20cm diameter, 0.1–100mTorr neutral fill, selectable gas and
Helmholtz-like adjustable bias are provisional exploration variables. Closed
flux or axisymmetric response does not establish 3-D stability/confinement.

---

# Continue from the wall-space correction and matched continuum gate

The full experiment-sandbox goal remains active. Read the newest GOAL_STATUS.md,
STRONG_WALL_DESIGN.md and STRONG_WALL_VERIFICATION.md. **436 tests passed in
232.03 s.** All new simulation, full replay and regression-test jobs are terminal.
No local or remote run needs continuation or restart.

## Current evidence

- Original rapid mode is localized near the outer wall with incompatible
  normal flow/Btheta/current traces. Independent physical reconstruction and
  radial/Fourier projections retain uncaptured norm and coarse background error.
- New explicit wall Galerkin curl Fw=P F, inactive wall-normal fluid momentum,
  admissible magnetic Krylov vectors and nonlinear Gauss stage coordinates.
  Strict new equations, full old raw-curl defect and zero reaction work are
  explicit. Original source classes and archived failures are preserved.
- Actual 344-coordinate contact growth 9203.373408/s matches the separately
  restricted original Jacobian. Three full matrices and all physical/projection
  arrays replay bitwise. Source ded3bac43979f5bde181c2eb6989d6af28b57f1d64bca09dd3f94413b7203a9a.
- Five new complete 100-microsecond contacts and six full waves; source
  d7a422b0f3e147e84e1acbb5b3db5c17b4eb162eab5d2d55845d829e18ebbc96.
  Wave time ratios 15.967079/15.991774; contact 15.883252. Contact energy drift
  decreases fourth order. Finest contact n error 9.97995e-7, pressures ~3.8e-6,
  E 5.00415e-4; unwanted B 1.00896e-13 T retained. All arrays/non-time metrics
  replay bitwise, with quadrature 16/20 audits.
- Independent original-continuum growing mode rate comparisons at mesh 1/2/4/8/16,
  quadrature 16/20 and amplitude 1e-4/1e-5: all 12 cases replay bitwise. Source
  c100c3eea05acbb767c65887bf692335ff88b788d884eef3839aaf7e5b0ec4fb.
  Finest full E rate comparison error 0.00648538, radial velocity 0.002272;
  radial-B/current rate errors remain 0.547444/0.480471. Do not claim adoption.
- 22 prior verified source closures and 32 previous-checkpoint archives unchanged.
  All 32 new source archives verified. Inspected scientific report PNG/PDF:
  outputs/strong_wall_report; checkpoint: outputs/strong_wall_diagnostics.
  Preserve user's modified outputs/frc_showcase/poster.png; no remote/publication.

## Immediate next action — frc-7mm.7 remains open

1. Read companions for continuum_contact_mode, continuum_wall_probe,
   strong_wall_geometry/plasma/induction, wall_induction_core and Gauss coordinates.
2. Refine radial and axial field meshes separately for the independent continuum
   mode. Current fixture couples default field resolution to fluid elements;
   generalize in a new explicitly frozen protocol, retaining analytic state,
   full physical rates/E/J, initial error and axis/wall/max/quadrature diagnostics.
   Slow magnetic response requires cancellation accuracy far beyond a B-field
   representation or fast-wave check. Do not suppress the residual mode.
3. Independently evolve small continuum perturbations with time/amplitude/mesh
   controls; compare physical fields/growth, not differing largest eigenvalues.
   Complete local charge/enthalpy/Poynting interface budgets including jumps.
4. Then continue heated moving exact vacuum, transport/ionization, actual coil/end
   boundaries, converged formation and experiment UI. Keep the complete goal and
   numerical acceptance beads open until their actual criteria pass.

This is substantial progress, not a blocked goal or a reduced objective. Older
sections below are historical; do not restart old terminated inaccurate jobs.

---

# Continue from implicit verification and exposed contact-growth gate

The full experiment-sandbox goal remains active, with twelve unfinished beads.
Read the newest GOAL_STATUS.md section, IMPLICIT_PLASMA_NUMERICS.md and
IMPLICIT_PLASMA_VERIFICATION.md. **432 tests passed in 227.47 s.** All jobs from
this checkpoint are terminal, including deliberately stopped inaccurate
contact candidates. No local or remote job needs continuation/restart.

## Current evidence

- Midpoint/Gauss/Radau advance the actual full physical-mass plasma equations
  with unchanged original induction bounds, fixed physical coordinate scales,
  strict nonlinear equations and atomic rejection. A fixed-background block
  preconditioner plus optional matrix-free Newton solves the equations without
  substituting the background into evolution. Dense coupling storage is only
  suitable for bounded verification, not large formation grids.
- Eleven full waves replay all arrays/non-time metrics bitwise; independent
  quadrature16/20 includes full E, both temperatures and physical energy.
  Time ratios: midpoint3.975482; Gauss15.967079/15.991775;
  Radau7.975201/7.994196. Radau has declared measured numerical temporal damping.
  Source feeb358468dfd0074e61effcca9af8dcb727748f113907d22b4e0e74e958c8ef
  (Gauss/midpoint),445623fdadd9f0825ae5441dc4c3fe26e92859e5c5d53564d5639949c9bb4240
  (Radau). Report: outputs/implicit_plasma_report/convergence.png and PDF.
- Three long Gauss contacts and one Radau contact developed large errors and
  were stopped, with terminal records and saved accepted states. Do not claim
  an accurate full100us transit. Initial chord-only failed pilots remain frozen;
  the captured ninth-step Gauss failure now has a passing Newton regression.
- A 360-coordinate stationary-contact spatial Jacobian has real growth
  854532.963/s at epsilon1e-5 versus854532.845/s at1e-4. The uniform rest control
  has maximum real4.28e-9/s. All matrices replay bitwise. Full source
  b82a7bf386a9d591e67699db71c9a5326e35c4bf709a30364b208645edd2449f.
  The uniform null-mode relative mixed error is ill-conditioned; separate
  nonzero independent rest-operator comparisons are <=4.403e-12.
- Paired species-pressure/current-entropy feedback is needed for rapid growth
  in saved-matrix ablations. Removing it is a diagnostic, not a solver fix.
- A separate continuum Fourier/Bessel operator imports no production spatial
  solver, verifies three uniform wave branches and has slower growing contact
  modes. First-mode axial refinement converges to8763.749254/s; higher sampled
  Bessel modes give about12,000/s. Eight matrices replay bitwise. Audit source
  b7dde0bbad1b307e752bc2a9b441e158dd7c790dd88728c684782667aee647a2.
  **The radial spaces differ: no matched-mode or numerical/physical cause has
  been established.** Spectral growth alone is not a laboratory stability claim.

## Immediate next action — frc-7mm.7

1. Read companions for contact_linear_probe, stationary_contact_continuum,
   stationary_contact_audit, regular_species, species_projection, regular_plasma
   and the relevant independent physical sampling code. Saved production modes
   and continuum matrices contain everything needed for the next comparison.
2. Reconstruct the production leading eigenvector as physical density, velocity,
   both entropy/pressure fields, B/J and local charge/face jumps. Compare its
   matched radial/axial subspace against the independently linearized parent
   equations. Consider an explicitly divergence-free continuum magnetic basis
   to reduce nullspace conditioning; the current full-field divergence residual
   is measured, not assumed zero. Independently evolve small perturbations and
   refine mesh/amplitude/time to distinguish physical growth from spatial error.
3. Fix only a demonstrated numerical defect. Do not suppress growing modes by
   dropping species feedback, fitting E/B/pressure, weakening residuals, adding
   artificial mass or repairing heat/energy. A stable uniform spectrum or fixed
   current species test does not certify the coupled heterogeneous state.
4. Repeat full contact and other unforced parent trajectories; close local
   enthalpy/Poynting interfaces. Then continue heated moving exact vacuum,
   transport/ionization, actual ends/coils, converged formation and experiment UI.
   Keep frc-7mm.4/.6/.7, frc-7mm and frc-goal open until their actual gates pass.

All previous source archives remain intact. No remote or publication; preserve
the user's modified outputs/frc_showcase/poster.png. Older sections below are
historical. This checkpoint made substantial new progress; it is not a blocked
goal and must not be resumed by resetting its scope or rerunning old remote jobs.

---

# Previous coupled moving-plasma checkpoint

The user made the full experimental-sandbox objective an active goal. Read the
new top section of GOAL_STATUS.md, REGULAR_PLASMA_NUMERICS.md and
REGULAR_PLASMA_VERIFICATION.md. The older isolated-species entry below is now
historical. The parent goal and eleven unfinished beads remain open.
**425 tests pass in 219.19 seconds.** Every job from this continuation is terminal;
there is no local or remote simulation to resume. Continue from saved evidence.

## Current evidence

- Actual fluid/current/two-species/magnetic RHS and full simultaneous RK4:
  25 frozen cases, complete bitwise replay, independent quadrature 16/20.
  Source 19f68fb08f949c4526bc1347d18652fb1eeb22d720d483a2074d6a71f7f40358.
- Seven full fast-wave periods; time-refinement ratios 15.999437/15.999863.
  Finest B/J signal errors 2.35089e-6/2.59599e-6; density/pressure 5.6764e-4.
  Contact has convergent original instantaneous rates only, not a full transit.
- Seven independent stable physical-energy audits, replayed at both quadratures;
  six slow-branch amplitude controls retain a weakest-pressure precision floor.
- The old full coordinate-matrix attempt's rejection was captured/replayed.
  **It is a localized near-axis toroidal magnetic perturbation, not the initial
  pressure-only hypothesis.** Inner physical and outer coefficient norms differ.
- Opt-in RefinedRegularPlasma tightens only the inner induction solve until all
  original unchanged outer residuals pass. Two full 360-coordinate matrices
  replay; independent linear errors 2.04814e-9/2.20696e-11. A complete coarse wave
  exactly matches the original. Source bc562813c4df25bb032e693b4e59adf0782aa9e92ee29f1233a36b34a1a2cd3f.
- Seventeen earlier calculation closures remain unchanged. All archives include
  complete source, arrays and metrics; do not overwrite or mutate frozen sources.

## Immediate next action

1. Read companions for regular_plasma, regular_plasma_step,
   refined_kernel_induction and refined_regular_plasma. Use the **opt-in refined
   class** for new work; original RegularPlasma remains the frozen baseline.
2. Complete a full exact FiniteMassContact transit on actual moving state with
   physical masses, both temperatures and full E. Derive a stiff time method or
   justify the explicit fast-mode cost; the contact transit is 100 microseconds,
   far longer than the 116.5 nanosecond fast-wave fixture. No rate clipping,
   artificial mass, field/pressure fit, heat or energy correction.
3. Preserve strict original residuals and atomic rejection across all stages;
   refine time/mesh and measure point/sampled positivity, initial/final physical
   norms and grid versus continuous conservation. frc-7mm.6 remains open until
   broader moving robustness/adoption is supported.
4. Close the original local enthalpy/Poynting budget with reconstructed face
   jumps. Integrated energy and nonlocal entropy graph work do not prove it.
   Inspect constant-potential gauge/EOS difference precision for very weak slow
   pressure signals; amplitude 0.001 controls retain a measured accuracy floor.
5. Only then continue heated moving **exact vacuum**, physical transport/
   ionization, actual coil/end boundaries, converged formation and experiment UI.
   Keep the original full goal open, including stability/confinement limitations.

There is no remote or publication. Preserve the user's modified
outputs/frc_showcase/poster.png. Existing study figures are verification reports,
not a proposed laboratory experiment. No old remote job requires continuation.

---

# Previous regular species checkpoint

The user resumed on 2026-10-04 and supplied a provisional lab envelope; research
is active and all ten previously unfinished beads remain open. **416 tests pass**.
Read the latest section of GOAL_STATUS.md, REGULAR_SPECIES_NUMERICS.md and
REGULAR_SPECIES_VERIFICATION.md. The species component is now implemented;
the earlier instruction below to begin that component is historical.

Frozen source `45a71e5b59cbcab6bb8556d531510d68da8389b16abcc82e1e61c0bc3237b78d`
has 24 bitwise-replayed cases, independent continuous audits, a local enthalpy
audit and eight full-transit fixed-current entropy trajectories. Actual moving
fluid/current evolution has not yet adopted it. Fifteen earlier calculation
sources are unchanged, and every job from this continuation is terminal.

## Immediate next action — frc-7mm.4 remains open

1. Read companions for regular_species, species_projection,
   regular_species_transport, regular_coordinate_fluid, regular_current,
   current_transport and kernel_induction before assembling a new coupled module.
2. Use the actual n,q,l,mz,Ce,Ci state and one common density derivative. Add
   species fluid_rate to bulk/current-stress rates; add **reduced_forcing** to
   induction, restoring **potential_electric** afterward for the full original
   electric/Ohm diagnostic. Do not numerically curl a large compatible gradient.
3. Verify independent original parent rates and complete unforced compressible
   wave/contact trajectories, retaining physical masses, local temperatures,
   axis/wall/max norms, time/mesh refinement and atomic failure. Derive the
   coupled time-step/implicit method for the actual fast electron/current modes;
   fixed-current RK4 stability does not select a safe coupled time step.
4. Close the local physical enthalpy/Poynting interface budget. The new graph
   transport is conservative and globally supported; it does not prove exact
   geometric face currents. The piecewise strong enthalpy divergence retains
   distributional jumps, and continuous energy remains separate from grid work.
5. Retain original fluid accuracy/spectral limitations and proceed to heated
   exact vacuum, physical closures, actual coil/end conditions and converged
   formation only after the coupled numerical gates pass. Keep frc-goal open.

The user's prospective 1 m quartz tube, 10–20 cm diameter, selectable gas,
0.1–100 mTorr neutral fill and adjustable Helmholtz-like bias are recorded in
RESEARCH_QUESTIONS.md. The verification fixture is not a selected lab experiment.
There is no Git remote. Preserve the user's modified outputs/frc_showcase/poster.png.

---

# Historical new-thread handoff

The user stopped this thread on 2026-10-04. The complete portable handoff is
[THREAD_HANDOFF.md](THREAD_HANDOFF.md), including all ten unfinished beads,
original scope, verified checkpoint, execution details and remaining gates.
No species implementation was started after research commit `7b3c8b2`.

The previous continuation note below is retained as historical context.

---

# Continue from the 2026-10-04 checkpoint

The user resumed on October 4. The overarching FRC goal remains open; this is
verified component progress, not a completed formation simulator. This repository
has no configured remote. Research and issue state are saved locally.

## Verified now

- Current advection, drift stress, stretching and self-drift: 16 frozen cases,
  bitwise replay, independent order-16/20 continuous-field audit. **408 full
  tests passed in 228.47 seconds.** Source `2f682eece854d54b3b0d33026709f8d34617cc4633d8b3875a0f3e2b8e0d9540`.
  The study CLI is complete; the missing records/audit modules from the prior
  pause have been implemented. frc-7mm.3 is closed. Read CURRENT_TRANSPORT_NUMERICS
  and CURRENT_TRANSPORT_VERIFICATION for exact scopes and measured limits.
- The old remote source-3f nr32 contact finished all 4116 steps and passed its
  original-platform final-step/physical audit. All 14 face-radial trajectories
  are complete and audited; frc-5b5 is closed. Evidence is in
  `outputs/face_radial_remote_completion/`. No restart occurred.
- Four species-pressure/entropy charge-compatibility probes, both quadrature
  orders replayed bitwise. Source c6a23e1a, `outputs/species_charge_probe/`.
  This diagnostic is separate from the 408-test solver baseline. It is not an
  implemented species transport operator.
- All 13 preceding calculation source closures are unchanged. Do not mutate
  their archived source definitions or rerun into existing immutable records.

## Next action: frc-7mm.4

Read `REGULAR_SPECIES_COUPLING_DESIGN.md` and the relevant module companions.
Derive the shared current/particle flux and both species entropy updates using
the actual regular n,q,l,mz,Ce,Ci state and physical h_e,h_i energy gradients.
The identity grad(p)/n=grad(gamma*h)-h*grad(s) exposes the required charge-flux
pairing. Compatible potential gradients annihilate current work, but their local
approximation error is nonzero and must converge; a global cancellation alone
is insufficient. Do not silently change density continuity or repair heat.

First implement and verify original local pressure and entropy rates with
nonuniform density, independently varying species temperatures and nonmatching
element traces. Measure axis/wall/max errors, species entropy, charge/particle
continuity, torque and physical versus grid energy. Then combine the verified
components and run unforced complete-parent wave/contact trajectories, including
enthalpy/Poynting boundaries and convergence. The lower-order continuous-energy
and earlier fluid order/spectral limitations remain recorded under frc-7mm.1.
frc-7mm.2 remains open for full coupling; frc-7mm.5 covers later adoption/formation.

## Later formation gates

Preserve exactly zero exterior plasma in heated moving-vacuum runs. Physical
transport/ionization, ends/coils and converged formation remain required.
Hardware stays abstract prescribed waveforms. Closed flux or axisymmetric
equilibrium does not establish 3-D stability or particle confinement.

## Machine and presentation

The authorized remote host is `m@192.168.0.202`; the completed old study is in
`/home/m/Code/frc-tool-comparison/face-radial-20260930`. The original PIDs were
absent on resumption and terminal evidence was audited/retrieved. No remaining
remote run is expected from that study; revalidate actual process identity
before any future action. Existing host authorization includes data transfer.
It was CPU NumPy/SciPy work, not a GPU-accelerated solver.

The six-panel convergence report is in
`outputs/current_transport_diagnostics/convergence.png` and `.pdf`.
The earlier interactive field presentation is preserved in
`outputs/frc_showcase_checkpoint/` (HTML, poster, MP4). It still depicts the
saved coupling case 07, with camera orbit only, and its historical status label
is from the previous checkpoint. Rebuild/serve via
`.venv/bin/python -m examples.showcase_checkpoint --serve --port 8779`.
The user's modified `outputs/frc_showcase/poster.png` remains outside our commits.
