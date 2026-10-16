# Active research: a credible experimental FRC simulation

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

## Current boundary: pulse/control physics and stronger prescribed field (2026-10-14)

Simulation-first continuation of3a471e54. Read MATERIAL_EXPERIMENT.md and
material_split_current.md. The complete bias-only metre-scale control63566 exits0:
final radius54.2617mm versus54.3501mm with the1mT opposing pulse. The pulse gives
field reversal and much greater heating, but almost no change in expansion.
Inspected actual Bz/flux/current/heat/radius figure:
outputs/material_tube_pulse_control_comparison/pulse_vs_bias.png (also PDF).

Third weak-pulse spatial run42299 exits0 at2000steps. Compared with6x12,8x16 has
final B/J/velocity differences10.05%/48.63%/0.3684%; local current remains unresolved.
Independent exterior/volume coil linkage discrepancy0.29949%. Actual-basis plot
and data: outputs/material_coil_spatial_r6_r8_q10. No additional test framework.

Initial central pressure157.936Pa corresponds to19.923mT magnetic pressure. A
20mT prescribed opposing-pulse exploration completed2000steps/3.5us (52842 exit0,
1197.04s), same assumed pre-ionized plasma. Minimum recorded midplane radius is
31.3557mm at2.10us, then53.1791mm at3.5us; maximum boundary radius68.6460mm is
larger than the control54.2617mm. Heating85.5213mJ, external work75.3002mJ,
meanTe/Ti7.12060/.439237eV. Final Bz-5.28416..+.929770mT with+.1mT source bias.
This is transient central compression plus greater outer expansion, not sustained
confinement.
Bulk particle-weighted RMS radius ends25.6110mm vs27.4959mm control (6.85% smaller),
both above18.7082mm initial. Thus the boundary maximum alone hides a modest bulk
effect; neither moment nor surface geometry establishes particle confinement. Drift/thermal ratio reaches116.93, exposing an invalidly established
edge closure; strong-run local/time accuracy remains unresolved. Work residual
maxima/E0 are2.7201e-5/4.9576e-5. No experimental/hardware recommendation.
Completed output: outputs/material_tube_bias_pulse_20mt_r4_z16_t2000.
Inspected comparison: outputs/material_tube_20mt_control_comparison.
Updated actual3-D/signed sections: outputs/material_experiment_20mt_display.

An opt-in split step handles the matched upwind pressure switch implicitly while
retaining all physical terms and directly measured external work. Two early
failed attempts and sources are preserved. Corrected candidate complete pulses:
COMPLETED14252: outputs/material_coil_split_r4_t2000_d,2000steps/7.121769us,
1066.17s. Actual time halving now reduces late sensitivity: total B/J differences
.21196%/.38724%, Jr/Jtheta/Jz3.7249%/.018683%/4.0554%, versus original RK4 weak
Jr/Jz15.28%/14.99%. Btheta still6.0173%; local accuracy is not fully established.
Max direct-work residual/E0 2.65160e-10, minima entropy within2.91e-13 of initial.
Inspected outputs/material_coil_split_time_comparison. No further new job pending.
COMPLETED45710: outputs/material_coil_split_r4_t1000_d,1000steps/7.121769us,
539.42s. Max direct-work residual/E0 1.06014e-9; minimum entropy deviation from
initial at most1.60e-13. Heat5.49400nJ, final peakJtheta9.91443A/m².
Actual method comparison with original1000-step RK4: final B/J differences
.9814%/1.8543%, weak Jr/Jz19.80%/18.11%; not a time-convergence certificate.
Inspected outputs/material_coil_split_vs_rk4_t1000. New step remains opt-in for physical cases beyond this weak pulse.
Both source d56d01a8aa762d4b58f241fb8abd3b5636278de928764c24497fdca1de1a0eeb.
Previous c runs exited1 at steps1656/808: a newly hit pressure bound repeatedly
released caused a zero-length solve step. The new retained working-set solve
resolves the actual captured quadratic to4.34e-19 without changing its equations.
The two complete weak pulses demonstrate improved time behavior; strong pulses
and negative pressure branches are not thereby validated. No unit tests.
LIVE89864: outputs/material_vector_free_r6_p5_q10_t2000, older full-vector run.
Poll exact handles; keep all live source closures unchanged.

Newest3-D viewer is localhost8790/20mt/ on server89915, browserTab8 marked
deliverable. Time/component/motion-label checks passed; no browser errors.
Scientific PNG/PDF and browser_view.jpg inspected. The old1mT view remains at8790/. Preserve the old views and user-modified poster. No remote.
The overarching goal and frc-vacuum.3 remain open. NEXT improve actual coil/transport/preparation simulations using these complete
results, including independently prescribed bias/formation coil groups. Use the
new step in concrete physical trajectories before broader adoption. Collect89864
when it exits; no repeat exhaustive grid/replay infrastructure.
Avoid expanding verification infrastructure. Keep the full goal open.

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

## Simulation-first continuation (2026-10-14)

The prior scientific turn made progress: F32 physical audit and full F8/F16/F32
comparison, with complete replays, finished and actual exits0 were collected.
Read COUPLED_FLUID_REFINEMENT_RESULTS.md. Physical energy-rate error decreases
with fluid resolution, but fluid response is not converged and rapid magnetic
response persists. No cause or cure is established by directional projections.

The user asks to prioritize the sim and end-to-end physical baselines rather
than further verification/reporting infrastructure. An actual F32/64,B32/64
100us/60-step negative1e-4 trajectory is now running through the existing
driver: parent session78724/PID33604, child33608. Paths:
outputs/trajectory_fluid32_campaign and outputs/trajectory_fluid32_t60_minus.
Poll the actual session; do not restart based on old status snapshots. Original
equations, bounds, duration and physical parameters remain unchanged. Storage
preflight reserved133413940 bytes per complete run with10.9GB free. No unit tests.

The full lab-design goal remains open. The immediate gate is accurate moving
plasma evolution; exact vacuum, transport/ionization/validity, physical coil/end
coupling and converged formation remain. The latest display still shows the
older accepted F8 step17 until a new physical trajectory is evaluated.

## Latest requested display and scientific boundary (2026-10-12)

Updated viewer: outputs/current_state_display_v3/index.html and inspected
current_state_3d.png. Signed Jr/Jtheta/Jz controls and quantitative meridional
Jz/Bz-departure slices accompany the unchanged full magnetic paths. Four real
saved frames still end at step17,28.3333us; step18 rejected, no reversal/FRC.
Complete builds replay byte-identically and all560 manifest inputs pass.
Browser interaction remains unverified because the Mac is locked. No unit
tests or physical solver runs for this presentation. Read CURRENT_STATE_DISPLAY_V3.md.

The interrupted scientific continuation has complete transferred directions,
F16/F32 signed-RHS probes/replays and the F16 physical audit/replay. These do
not provide a later evolved state or demonstrate a growth cure. F32 physical
audit and the prepared comparison are not executed. Read
COUPLED_FLUID_REFINEMENT_PROGRESS.md for exact evidence, expired-session
qualification and next command. Do not restart completed calculations.
Full FRC goal and all acceptance gates remain open. Older next-action blocks
below are historical and superseded by this boundary.

## Latest continuation: entropy budget and fluid-space resolution (2026-10-12)

Verified progress from 1146f0b. Read SPECIES_INTERFACE_VERIFICATION.md and the
inspected outputs/species_budget_report/verification.png. The whole-vector
independent weak entropy identity matches within 1.78627e-12 relative, including
all one-sided faces and projection adjoints. Its original production comparison,
analytic control, both species, both directions and q12/16/20 fully replay.
Large interface/transfer terms cancel; they cannot individually explain growth.

Independent physical mass diagnostics leave most of the 55% strong-source error
as an unresolved representation remainder. With original B32/64 current fixed,
F8/F16/F32 refinement reduces original electron source error
0.5500115 → 0.0897494 → 0.000292433; ion behavior is similar. Exact-mass vectors
remain diagnostics only. This is scalar-source convergence, not a demonstrated
repair or coupled growth/trajectory convergence. The complete 128-row actual
q16/20 field comparison has maximum full-rate-relative difference 1.93384e-12.

All eight calculation/report workers are terminal; actual exits0 collected.
The 9+4+12+1 raw files, complete metrics and report PNG replay exactly. Six
guarded jobs peak at 696647680 bytes and retain at least 5939200000 free bytes;
report stage observations are separate. No unit tests or physical solver change.
Both display versions, frozen sources, original states and user poster remain.
Full goal and all acceptance beads remain open.

Final complete 200-source preservation record:
outputs/species_budget_diagnostics/checkpoint_v2.json.

NEXT execute COUPLED_FLUID_REFINEMENT_DESIGN.md: exactly subdivide both full
diagnostic directions onto F16/64,F32/64, verify all physical transfers, and
evaluate the unchanged complete signed RHS with original residual bounds and
amplitude controls. Then choose full signed/time/arithmetic trajectories and
local budgets before vacuum/transport/coils/ends/formation gates. No old physical
worker should be restarted. Earlier status sections below remain historical.

## Requested display update (2026-10-12, after d4cff11)

The current 3D viewer and inspected PNG are now in
outputs/current_state_display_v2. Read CURRENT_STATE_DISPLAY.md. The final
accepted state remains step 17 at 28.3333 µs; later audits have not advanced it.
New camera presets, keyboard controls and latest audit status accompany the
unchanged saved B/J data. Complete builds and their fresh-process replays are
byte-identical; no unit tests or solver runs. Browser inspection remains
unavailable while the Mac is locked. The first clipped-camera PNG and all older
sources/posters remain preserved; the delivered v2 has unclipped geometry.
The full research goal and all acceptance beads remain open. Next scientific
action remains SPECIES_INTERFACE_DIAGNOSTIC_DESIGN.md below.

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

## Refined amplitude and residual controls complete (2026-10-12)

All workers are terminal and collected. No unit tests were added or run.
This continuation from cd66ec5 made verified progress; the full FRC goal and
all12 acceptance beads remain open. Read REFINED_CONTROL_DESIGN.md and
REFINED_CONTROL_VERIFICATION.md. Final checkpoint:
outputs/refined_control_diagnostics/checkpoint.json.

On each independently refined field grid, the full epsilon1e-5 and1e-3 signed
workflows pass the unchanged1e-10 original checks and fresh complete replays.
Both epsilon1e-4 / tolerance1e-11 controls strictly reject on minus after all
five outer attempts, then reproduce their complete rejection protocols.
Final dual-Faraday residuals are6.14598e-11 and8.29280e-11. Each rejected case's
five saved input arrays match the accepted baseline in shape, dtype and bytes;
this does not establish post-call in-memory immutability. No failed rate is used.

Four component, two spectrum and two actual magnetic-rate curl audits fully
replay and reproduce original physical/reference/weight fingerprints. On
B256/320, Jr error rises1.76199% to4.17506% at the smaller amplitude; independent
strong curl gives1.73644% to3.19274%. On B320/256, Jr rises5.90438% to6.70354%.
Actual small/baseline Jr differences are3.83376% and3.19588%; total-J norms hide
them. The larger amplitude slightly worsens radial/axial current, velocity and
azimuthal-field errors. No amplitude-convergence or selected-waveform claim.

Above cell mode|k|64, radial-current error grows9.9952x and9.7616x as epsilon
shrinks10x. That band holds93.7982% and87.2576% of squared amplitude difference.
This is consistent with numerical noise amplified by differencing, without
identifying its source. Significant lower-band error remains. Cell modes fold
aliases; strong curl excludes face distributions. No state or mode is filtered.

All240 independent component/spectrum integrals and48 component/curl integrals
match bitwise. Volume component q16/20 norms agree within3.440e-15, while
band-norm differences reach4.968e-5 and remain separately reported. Parseval
and curl decomposition defects stay below3.150e-16. Instantaneous exchanged-
power closure is within3.983e-8; local/time-integrated budgets are still open.

All12 physical protocols satisfy8GiB/256MiB guards including final extrema;
maximum peak7.574GiB, minimum free disk605.469MiB. Audit observations are1.464GiB
peak and294.531MiB minimum at stage boundaries, not continuous guards. No
interruption occurred. All byte-identical shared replay paths remain immutable.

Fourteen replay groups and14 raw-file comparisons are complete. Four new
modules have companions;58 new archives,446 preceding archives,164 preceding
sources,228 earlier evidence files and22 older active closures verify. Final
168-file source SHA:5deb5dc6c4a3229d60d13d33b2a182716aadd61fcd1fcfa1fa5264ef6dff4966.
Scientific PNG/PDF: outputs/refined_control_report; PNG visually inspected.
Original physics/defaults/time integrators and user poster are preserved.

NEXT resume full signed radial1/N64, F8/64,B32/64 100us trajectories with60/120
steps, independent physical/conservation audits and complete stage/state replays,
against the preserved F4/64,B16/64 controls. Budget resources and retain any
failure or fast-growing response. This is separate from the radial16 rate study.
Keep amplitude noise/strict residual arithmetic and joint B320/320 feasibility
open. Local budgets, heated exact vacuum, transport/ionization/validity, real
prescribed coils/ends, converged formation and UI remain required. Lab inputs
stay provisional; no3-D stability or particle-confinement claim. No remote or
publication. Earlier completed sections below retain their historical scope.

## Separate field refinements verified with lower workspace memory (2026-10-12)

All workers are terminal and collected. No unit tests were added or run.
The full original FRC sandbox goal and all 12 acceptance beads remain open.
Read SCHEDULED_WORKSPACE_DESIGN.md and SCHEDULED_WORKSPACE_VERIFICATION.md.

Complete coarse/fine signed workflows and fresh replays establish bitwise
original/new fidelity for scheduling and componentwise drift-stress evaluation.
The unchanged stress adjoint receives each original component expression on
request, retaining multiplication association and all physical checks. Fine
peak memory drops from scheduling's 6.481/6.476 GiB to 5.342/5.516 GiB.
Default solvers/time integrators and old calculation sources remain unchanged.

At fixed F64/256, radial16/N128 reference, epsilon 1e-4 and tolerance 1e-10,
independent axial B256/320 and radial B320/256 signed calculations, physical
audits and full fresh replays pass. Axial refinement reduces Jr error from
5.96939% through 3.05317% at axial288 to 1.76199%; independent strong curl gives
1.73644%. Br error falls to 0.0266691%. Radial refinement instead reduces Jz
error from 1.18700% to 0.734849%, and radial-velocity error from 0.183143% to
0.0983390%; Jr stays near 5.9%. Separate directions address different errors.
The modest increments do not establish asymptotic or amplitude/time convergence.

New maximum residuals 6.82535e-11/8.96429e-11 obey 1e-10. Instantaneous work
closure is within 4.343e-9 of exchanged powers; local/time budgets remain open.
New physical q16/20 rate norms differ by at most 1.39e-15, component norms by
2.64e-16; curl decompositions close within 2.30e-16. Strong curl excludes face
distributions. No field, mode, state or strict acceptance bound was changed.

Axial primary/successful replay peak 6.901/7.162 GiB; radial 7.161/7.037 GiB.
An earlier axial replay was disk-guard interrupted and collected as exit130;
its full partial evidence remains. A fresh standalone replay passes. All
successful final resources satisfy 8 GiB RSS / 256 MiB disk guards. Verified
byte-identical completed NPZs share immutable storage; every path/hash survives.

Twelve replay groups: six full signed simulations, four implementation
comparisons, two independent curl audits; eight raw-file comparisons. Eight new
modules have companions; 49 new archives, 397 preceding archives and 22 older
active closures verify. Final 164-file source SHA:
a7022b5b4d1a8f5b598c431e71f6865d51c178ca61eed43a33ef230cfaeeabf5.
Checkpoint: outputs/scheduled_workspace_diagnostics/checkpoint.json.
Scientific PNG/PDF: outputs/directional_workspace_report; PNG visually checked.

NEXT add amplitude/residual controls on the two refined grids, then resource-
budget combined/further mesh refinement with component/actual-curl norms.
Full signed trajectories/local budgets, heated moving exact vacuum, physical
transport/ionization, prescribed coils/ends, formation and UI remain open.
User poster and all evidence preserved; no remote/publication. Lab inputs stay
provisional; these filled-periodic fixtures prove neither 3-D stability nor
particle confinement. Older sections below are historical.

## Axial refinement halves radial-current error (2026-10-12)

All workers are terminal and collected. No unit tests were added or run. Read
WORKSPACE_MEMORY_DESIGN.md and WORKSPACE_MEMORY_VERIFICATION.md. The full
original FRC sandbox goal and all12 acceptance beads remain open.

Two opt-in workspace changes preserve complete original coarse/fine physical
workflows bitwise: release expired quadrature temporaries, then detach completed
density models from cached geometry templates and collect unreachable inertia
cycles after each full rate. All four original/new comparisons and their replays
pass, and all complete candidate runs have fresh all-array replays. Default
geometry, old sources and time integrators remain unchanged.

The lifetime-only fine run/replay peak6.873/7.430GiB; cache detachment peaks
6.446/7.056GiB. This supports a modest guarded field-axial increment256→288
at fixed fluid64/256 and radial field256, independent radial16,N128 reference,
epsilon1e-4 and original tolerance1e-10. The new signed run and full fresh
replay both pass, at7.268/7.408GiB; no guard triggers. Doubled/radial-grid
feasibility remains unproved. Peak observations are not memory bounds.

Independent q16/20 physical component and magnetic-curl audits reproduce all
original fingerprints and replay completely. Jr error drops5.96939%→3.05317%;
strong-curl Jr error5.94247%→3.03377%. Weak/strong difference slightly grows
0.110758%→0.120095%, so improvement also lies in the actual magnetic response.
Br error0.0979624%→0.0491827%; Btheta0.547509%→0.363014%. Jz remains1.11402%
and the full-current vector error0.0817008% still hides weaker components.
Two grids establish neither convergence order nor full time/amplitude accuracy.

Worst refined original residual6.63361e-11 remains below1e-10. Largest
instantaneous power defect3.111e-9 of exchanged powers; no local/time budget
claim. Quadrature component norms agree within4.76e-16 absolute, physical rate
norms within1.30e-15; curl error decomposition closes4.23e-16 relatively.
Eleven replay groups comprise five complete signed calculation replays, four
whole-workflow comparison replays and two independent curl-audit replays.
Nine new modules have companions;42 new archives,355 preceding archives and22
older source closures verify. Final156-file source SHA:
7e89bbaac2a377a3abe1d2a6fab0b0641381dbd777ae04a903a5c86e8c32e918.
Checkpoint: outputs/workspace_memory_diagnostics/checkpoint.json.
Scientific PNG/PDF: outputs/workspace_refinement_report_layout; PNG inspected.

NEXT reschedule independent gradient construction after stress evaluation to
reduce peak live tensors, retaining every expression and all original checks;
verify full original/new signed workflows and replays before larger axial/radial
refinement. Add refined-grid amplitude/residual and local-budget controls.
Full trajectories, heated moving exact vacuum, transport/ionization, real
prescribed coils/ends, converged formation and experiment UI remain open.
All scientific evidence/partial attempts and the user's poster are preserved;
byte-identical completed replays share immutable storage with manifests. Lab
inputs remain provisional. No3-D stability/confinement claim or publication.
Older sections below are historical.

## Separable mass inverse verified; field refinement awaits lower peak memory (2026-10-12)

All workers are terminal and collected. No unit tests were added or run. Read
FIELD_COMPONENT_REFINEMENT_DESIGN.md and FIELD_COMPONENT_REFINEMENT_VERIFICATION.md.
The full original FRC sandbox goal and all 12 acceptance beads remain open.

A new component audit of existing F32/256, B256/256 full signed data replays
bitwise. Fluid refinement to F64 on the same field grid reduces radial current
error 11.6635% to 5.96939%, axial current 36.1204% to 1.18700%, and azimuthal
current 0.135383% to 0.0190466%. Strong-curl comparisons reproduce saved fields.

The unchanged B256/512 initialization was resource-interrupted before any
signed input or derivative existed. Own peak RSS reached 11.230 GiB; disk guard
crossed 256 MiB. Approved termination was collected as exit 143. Partial records
remain preserved; this is neither numerical rejection nor an accuracy result.
Source storage estimates exclude large LU fill and workspaces and are not peak
bounds. No system-wide disk cause or exact slow function is established.

An opt-in exact tensor-product mass inverse replaces global 2-D mass LU with
separate radial/axial solves, retaining original assembled physical equations,
all modes, actual masses/density, quadrature, boundaries and strict checks.
Complete coarse original/new signed workflows capture identical actual operator
fingerprints and replay their archived pilots. Full fine signed F64/256,
B256/256 runs and fresh all-array replays pass. Worst residual 6.4967e-11 under
1e-10; instantaneous power defect <=8.75e-9 of exchanged powers. Independent
component audits reproduce original fingerprints and replay completely.
Fine old/new radial-current physical difference is 1.04097e-5 of reference,
while its error stays 5.96939%. This is not a spatial-accuracy repair.

Fine baseline/replay peaks are 7.420/7.150 GiB. A complete read-only resource
observation replay peaks at 7.432 GiB with identical scientific results. First
current transport raises cumulative peak 4.829 to 6.558 GiB; the second signed
electromagnetic evaluation raises it 7.061 to 7.432 GiB. These nested high-water
intervals are not exclusive/live allocation measurements. No new guard fires.
Default geometry and time integrators remain unchanged. Larger-grid feasibility
is not established. Every source/data/partial result and user poster survives;
byte-identical completed replay NPZs share storage with explicit manifests.

Seven new modules have companions. All 37 new source archives, 318 preceding
archives and 22 older source closures verify. Final 147-file source SHA:
8ae692975fec4c2fac1345f723b13bb428c1bb1affa2c72557f527d476c8b332.
Full checkpoint: outputs/field_component_diagnostics/checkpoint.json.

NEXT reduce measured temporary quadrature storage without changing equations,
then full same-grid physical/work comparisons and complete replays before
retrying independent field-axial/radial refinement. Keep individual current
norms and actual Bt derivatives. Complete trajectories/local budgets, heated
moving exact vacuum, transport/ionization, real prescribed coils/ends,
converged formation and experimental UI remain open. Lab inputs are provisional;
no formation, 3-D stability or confinement claim. Local only, no publication.
Older sections below are historical.

## Curl arithmetic separated from persistent component error (2026-10-12)

All workers are terminal and collected. No unit tests were added or run.
Read CURL_EVALUATION_DESIGN.md and CURL_EVALUATION_VERIFICATION.md. The prior
checkpoint 398f9a3 made verified progress; this continuation preserves it.
The full experimental FRC goal and all 12 acceptance beads remain open.

Every original residual/lift vector is reproduced bitwise on the captured
strict rejection. Factoring the projected curl alone does not fix it. Exact
raw products/sums plus addition-low accounting reduce the distributivity gap
from 4.55e-11 to 2.97e-15 of RHS. The actual saved iterate still has dual residual
5.62e-11 against 1e-11. No solver, equation, bound or accepted state changed.
The full observation workflow and fresh reconstruction replay match bitwise.

Complete physical cell-mode audits at epsilon 1e-3/1e-4/1e-5 reproduce original
field/reference fingerprints and fully replay. All 120 previous component
integrals match exactly; Parseval closes to 2.03e-16. Radial-current errors are
5.9751%/5.9694%/6.6721%. Above-reference-band norm errors are 0.02750%/0.27966%/
2.78699%; dimensionless high-band norm times epsilon is nearly 2.75-2.80e-7.
This is consistent with a numerical noise floor amplified by differencing,
while roughly 6% in-band error persists. The |k|20 error hump is shared by all
amplitudes. No field or state is filtered; cell modes fold continuous aliases.

Independent differentiation of the actual magnetic rate yields 5.9425% radial
current error at epsilon 1e-4, close to Jt's 5.9694%. Weak/strong difference is
0.1108% of the reference. Signed cross terms close the error decomposition to
2.16e-16; strong curl excludes distributional face jumps. This directs work
toward magnetic-response accuracy, without declaring the weak Ampere law wrong.

The full signed epsilon 1e-3 F64/256,B256/256 rate calculation and independent
q16/20 audits complete and fully replay, including every coefficient and face
array. Worst original residual 1.91046e-11 remains below 1e-10. Larger amplitude
slightly increases velocity/toroidal-field error; it is not a convergence fix.
Five verified groups comprise one full signed rate workflow and four diagnostics,
not time trajectories or formation experiments. Six new modules have companions;
18 new archives, 300 preceding archives and 22 older source closures verified.
Final 140-file source SHA b0ee5f17e5baf0fa4fc564d3fe40437173b05479d16210e2ef040f0a1d0ce922.
Full checkpoint: outputs/curl_evaluation_diagnostics/checkpoint.json.
Figure: outputs/curl_spectrum_report_layout/verification.png, visually inspected.

NEXT budget separate axial/radial field refinement on the same independent
radial16,N128 reference, comparing individual current components and actual
magnetic-rate derivatives. Investigate accurate raw curl/electric representation
with independent original equations, adjoints/work and full strict replays;
current diagnostics are not a repair. Then complete full signed trajectories,
local budgets, heated moving exact vacuum, transport/ionization, prescribed
physical coils/ends, converged formation and the experimental sandbox UI.
Provisional lab inputs remain exploratory. No 3-D stability or confinement
claim. User poster unchanged; no remote/publication. Older sections follow.

## Component errors exposed and strict residual paths isolated (2026-10-12)

The full FRC goal and all 12 acceptance beads remain open. No unit tests were
added or run. All workers are terminal and collected. Read
SPLIT_RESIDUAL_PATH_DESIGN.md and SPLIT_RESIDUAL_PATH_VERIFICATION.md.
The previous turn made verified progress; this continuation preserves it.

A new complete physical audit reproduces every archived rate/E/current/reference
fingerprint and fully replays. At fixed radial16, N128 and F64/256, B256/256,
radial current-rate errors are 5.9694% / 6.6721% at amplitudes 1e-4 / 1e-5.
Axial-current errors are 1.1870% / 1.1936%; azimuthal-current error is 0.01905%.
The previous 0.0879% total-current error remains correct but hides these weaker
components. Actual radial-current amplitude difference is 3.0148% of its own
reference, with 88.674% of squared difference in the interior, 11.115% in the
outer 5% radius and 0.210% in the inner 5%. Radial-current error is about 6%
across all bands. By contrast, 58.08% of axial-current squared error and 72.05%
of radial-velocity squared error lie in the wall band. Error and sensitivity
fractions remain distinct. Whole-domain q16/20 norm changes are <=3.44e-15;
small masked-band quadrature differences are explicitly retained.

An observation-only exact-input workflow reproduces the archived tighter-bound
first-attempt rejection. The actual opposing residual-vector gap is 6.87234e-11
of the RHS after three corrections. Electric curl distributivity accounts for
that vector to a remainder of 0.4245%; H/T association and inertia-path gaps
are much smaller. All seven observed contributions leave only 5.41e-17 of RHS
unexplained. An independent Decimal audit verifies all 786432 electric-addition
coefficients exactly, but addition-loss propagation alone leaves 66.21% of the
gap. Curl arithmetic must therefore remain separate from electric addition.
No solver arithmetic, physical equation, bound or accepted state was changed.

The first in-process observation replay was stopped by its disk guard. Original
and partial outputs are preserved. A fresh standalone replay completed and
matches every array, iterate and non-runtime metric; its guard did not trigger.
The independent arithmetic audit and component audit also fully replay. Three
verification groups and 12 new source archives are checked; 288 preceding
archives and 22 earlier closures remain unchanged. Six new modules have
companions. Final 134-file source union:
150dc2f3e7614eed5da58d314b6ff3e8cf804b641c1f9fec9f91b738e4712745.
Full checkpoint: outputs/split_residual_diagnostics/checkpoint.json.
Scientific PNG/PDF: outputs/split_precision_report; PNG visually inspected.

NEXT compare assembled projected curl with curl followed by projection on the
exact captured vectors, without relaxing checks or claiming a repair before
full physical verification. Separately establish component-wise current-rate
spatial/amplitude convergence; their relationship to this arithmetic failure
is not yet established. Then resume full signed trajectory convergence and
local budgets, heated moving exact vacuum, transport/ionization, physical
ends/coils, converged formation and the interactive experimental sandbox.
The lab envelope remains provisional. These fixtures establish neither
3-D stability nor confinement. User poster unchanged; no remote/publication.
Older sections below are historical.


## Combined refinement, full replay and resolved amplitude metrics (2026-10-10)

The full FRC goal and all 12 acceptance beads remain open. User policy:
**end-to-end workflows only; do not create or run unit tests.** No unit tests
were added or run. All workers are terminal and collected. Read
JOINT_E2E_DESIGN.md and JOINT_E2E_VERIFICATION.md. Older counts below are historical.

At the fixed radial16, N128 independent reference, F32/256, B256/256 gives
radial velocity / radial B / toroidal B / current-rate physical errors
0.0308854 / 0.000982554 / 0.0494785 / 0.0259173. F64/256 on the same field grid
improves these to 0.00183143 / 0.000979624 / 0.00547509 / 0.000878993.
Both full signed cases and fresh reruns match every coefficient, physical
fingerprint and non-runtime metric. Original strict bounds remain; worst
residual is 6.4614e-11 under 1e-10. Quadrature 16/20 differences are at most
5.00e-15 absolute; a 10.44% relative axis-density outlier is an error near 4e-14.

The finest epsilon 1e-5 case and full replay also pass. Current/velocity errors
are 0.000885440 / 0.00183390. Centered current coefficients differ by 0.00947096
in their Euclidean norm. A new independent physical reconstruction reproduces
all original current/E fingerprints and fully replays. The actual vector
current difference is 0.000107450 relative to the continuum reference, but the
small radial component changes 0.0301482 relative to its own reference; axial
current changes 0.00126286. Full E changes 7.14889e-8. The total current norm
therefore hides a 3.01% radial-component sensitivity. Quadrature controls agree.

The tighter 1e-11 bound rejects the identical minus input. Its complete rerun
reproduces every failure message and saved array. Dual Faraday residual stays
near 5.89e-11 after all original budgets; no rate is accepted for that case.
The two evaluated residual paths differ by at least 2.499e-11 at the recorded
final iterates. This is a lower bound, not an identified cause or proof that
another iterate cannot pass. Tighter-residual convergence remains open.

Separately, both radial1, N64 signed 100 us / 60-step trajectories and full
reruns match every accepted stage/state: 62 NPZ files per sign. Independent
physical, conservation, old/new and signal-energy audits fully replay. Final
current/E errors are 0.00274791 / 0.00867839; old/new final fields agree exactly.
Energy drift is at most 3.636e-7 of initial physical perturbation energy, which
is a normalization scale rather than an asserted conserved free energy.

All earlier Python sources remain unchanged. Opt-in streamed audits, reused
axial bases, omission of an immediately discarded LU and immutable geometry
reuse preserve the recorded physical calculations. Complete original/new
comparisons and adopted-variant replays verify these changes. All interrupted
parent batches and partial replays remain preserved. The 20 verified groups
include full simulations, explicit rejection replay, physical postprocessors
and completed-case selection; they are not 20 distinct successful experiments.
All 71 new source archives, 217 preceding archives and 22 earlier source
closures are verified. Final source union: 128 files, SHA
8a97cac82e67e4aa62eeea9f5c80edfd1c982529a158edd6fac0512d9f56e457.
Full inventory: outputs/radial16_joint_diagnostics/checkpoint.json.
The scientific PNG/PDF is outputs/joint_e2e_report; the PNG is visually checked.

NEXT capture the exact tighter-bound split-inverse residual vectors and locate
the radial-current amplitude sensitivity. Then continue full signed trajectory
mesh/time/residual convergence and local budgets. The finest trajectory's
retained dense modal arrays alone require an estimated 15.414 GiB, excluding
field operators and workspaces; budget this separately from rate calculations.
A smaller F8/64, B32/64 full-transit comparison is specified in NEXT_SESSION.

Heated moving exact vacuum, transport/ionization, actual ends/coils, converged
formation and experiment UI remain open. Lab dimensions, gas, neutral fill and
bias remain provisional. These filled-periodic fixtures establish neither
three-dimensional stability nor confinement. The user's poster is unchanged.
Research and issue state remain local; no remote or publication. Older sections
below are historical.

## Directional error separated; strict arithmetic recovered (2026-10-10)

The full experimental FRC sandbox goal remains open. **470 tests pass in233.11s.**
Read DIRECTIONAL_TERM_DESIGN.md and DIRECTIONAL_TERM_VERIFICATION.md. All new
workers are terminal and collected:13 complete bitwise verification groups,
42 verified source archives. No worker needs resuming. The23 preceding radial16,
59 physical-interface,61 matched and32 wall archives are byte-unchanged;22
older source closures still match. No previously tracked Python source changed.

At the fixed radial16,N128 reference, independent axial refinement reduces
Br-rate error .051938 to .0023194. Independent radial-field refinement reduces
vr/J-rate errors .1088365/.0875688 to .0309084/.0337504; Btheta remains .2793811
and Br .051625. Fluid-only doubling improves radial pressure error31.53-fold,
but leaves the dominating magnetic-force error. Physical accuracy remains open.

Original parent-term audits expose Hall/pressure toroidal induction norms about
24000 times their small net response. All signed original fluid terms reproduce
the full archived rate bitwise. Fixed physical outer-wall-band momentum-error
fractions fall from89.11% to2.263% under field refinement, but rise to96.60% with
fluid-only refinement. Driving-only audits preserve the old failed full-solve
status and never invent an inertia or accept a magnetic/current derivative.

The exact radial256 old solver failure remains preserved and fully replayed.
Observation rules out the initial HT-association hypothesis. Mass-inverse
floating-point distributivity explains the dominant dual residual gap, retaining
7.06% unexplained vector remainder. A new opt-in split-inverse arithmetic computes
E=M^-1f+M^-1(IJt) while preserving the original equation, masses, actual density,
boundary space, operator, every residual definition/bound and all budgets.
The exact captured full input now passes at5.27777e-11 under1e-10 and replays
bitwise. Independent full power defect is3.99063e-9 of exchanged powers
(absolute2.50e-19W); full Faraday rate-relative2.722e-7 remains separately visible.
Analytic full cylindrical B/J/E/work and atomic-failure controls pass.

Three new full signed pairs, their full reruns, amplitude/baseline/quadrature
comparison and all term diagnostics are terminal and bitwise verified. The
same-grid original fluid rates remain bitwise identical; magnetic/current
coefficients change only4.39e-11/5.86e-10 relatively under new arithmetic.
Amplitude1e-4/1e-5 physical errors nearly agree, while their current coefficient
difference .00108516 is explicitly retained. The frozen full trajectory solver
has not adopted the new arithmetic. No diffusion, filter, atmosphere or repair
was introduced.

NEXT compare F32/256,B256/256 with F64/256,B256/256 using the unchanged independent
radial16,N128 reference, fresh archives, both signs and amplitude/quadrature/
strict-residual controls. Then establish complete signed trajectory mesh/time/
residual convergence and local face budgets. Moving heated exact vacuum,
transport/ionization, actual ends/coils, formation and experiment UI remain open.
The user's lab envelope is provisional; this filled-periodic verification study
is not a selected experiment and proves neither3-D stability nor confinement.

Inspected PNG/PDF: outputs/radial16_directional_report_v2. Full inventory:
outputs/radial16_directional_diagnostics/checkpoint.json. Preserve the user's
poster. Research/issue state are local; no remote/publication. Older sections
below are historical.

## Fixed radial wavelength and recovered strict solve (2026-10-09)

The full laboratory FRC sandbox goal remains open. **462 tests pass in 230.35 s.**
Read RADIAL_MODE_DESIGN.md and RADIAL_MODE_VERIFICATION.md. A new opt-in
balanced dual-residual induction solver changes only iterative refinement.
Every original equation, physical mass, boundary, residual bound and refinement
budget remains; no old calculation source or trajectory solver was changed.

All jobs from this continuation are terminal and collected. Seven full bitwise
verification groups and23 new source archives are verified. The59 preceding,
61 older matched and32 older wall archives are byte-unchanged, and22 earlier
verified source closures still match. No worker needs resuming. The full
checkpoint is outputs/radial16_diagnostics/checkpoint.json.

The independent radial16 reference now has N64/128/256 axial controls. Its
N64 growth was nearly converged but off-grid curl/current disagreed by .325%.
The resolved N128 reference agrees with N256 physical/J/E shapes within 6.9e-9
and passes independent finite-difference curl controls. All mesh comparisons
hold this one physical wavelength and reference fixed.

Separate fluid/field refinements reduce full current-rate error from 5.10160
at F16/128,B64/128 to .0875688 at F32/128,B128/128. Radial velocity improves
from 3.86912 to .108837; Br remains .051938 and Btheta .383069. Initial
representation and exact one-sided face jumps improve much faster. The physical
accuracy gate remains open; this is an instantaneous bounded periodic fixture.

The old finest solve rejected the exact saved input at operator residual
1.39069e-10 versus the unchanged 1e-10 limit. Balanced correction proposals
make all original residuals pass, worst 7.72779e-11, with bitwise exact replay.
Independent analytic polynomial B/J/E and work controls pass on three grids;
iteration exhaustion is atomic. The old successful baseline remains bitwise
identical. Amplitude1e-4/1e-5 errors nearly agree, with explicitly retained
coefficient-vector differences. Volume quadrature16/20 norms differ at most
3.05e-11 relatively; tiny axis errors and their absolute scales are documented.

Next separate axial refinement F32/256,B128/256 from radial-field refinement
F32/128,B256/128, then fluid64. Localize the remaining pressure/current/work
and wall/interface rate defect before repeating complete signed trajectory
mesh/time/residual convergence. Preserve moving exact vacuum, physical transport,
coil/end boundaries, formation and experimental UI as unfinished full-goal work.
The user lab envelope remains provisional; no lab experiment is selected.

The inspected scientific figure is outputs/radial16_report/verification.png
and PDF. The final source/replay/worker inventory is recorded in
outputs/radial16_diagnostics/checkpoint.json. Older sections below are historical.


## Local budgets verified; radial response accuracy fails (2026-10-05)

The full experimental sandbox goal remains open. **457 tests pass in254.80s.**
Read PHYSICAL_INTERFACE_AUDIT_DESIGN.md and PHYSICAL_INTERFACE_VERIFICATION.md.
No production evolution source, physical mass, strict tolerance or accepted
state changed. Every new diagnostic and physical audit replays bitwise.
Four new100us signed trajectories and all four full evolution reruns are
complete. **Every accepted state/stage and non-runtime metric replays bitwise;
all jobs are terminal.** Earlier sections are historical.

Eight independent original-equation cell/one-sided-face audits show energy
defect.604850 -> .0453847 -> .00281026 -> .000173539 with refinement; current,
particles and both species entropy also improve. All raw face/cell/jump arrays,
axis/wall norms and quadrature16/20 controls replay. Small net energy power does
not certify local balance. New thermal-Hessian tests give a physical initial
perturbation energy scale; actual total-energy drift is<=6.42e-7 of that scale
over five signed pairs. The scale is not a separate conserved free energy.

**Full growing-response mesh/time convergence still fails.** At F4/64,B16/64,
the intended mode projection is more accurate, but J/Bz errors rise to
.0431611/.00496073. Coarse30/60/120 steps do not show fourth-order convergence.
Saved-history extraction gives unwanted current growth71682/s; direct actual
stationary spatial RHS gives71426/s with amplitude controls. Independent radial
projection locates it near Bessel modes15/16, on16 radial field cells. Strong
volume curl(B)/mu0 differs from weak J by.599823; interface distributions and
uncaptured radial norm remain explicit. Exact fluid-only subdivision preserves
fields but strongly changes momentum rates; it is not a demonstrated repair.
Representative uniform spectra also expose unresolved fast temporal phase,
which does not explain away the directly measured spatial contribution.

NEXT: independently matched radial/field/fluid and fixed-physical-wavenumber
controls, then complete signed time/mesh/residual convergence. Keep frc-7mm.7/.6/
.4 and the full goal open. Heated moving exact vacuum, transport/ionization,
physical coil/end boundaries, converged formation and experiment UI remain.
The lab envelope remains provisional. Inspected six-panel PNG/PDF:
outputs/physical_interface_report. All59 new source archives are verified;
61 preceding archives,32 older archives and22 earlier source closures remain
unchanged. Full checkpoint: outputs/physical_interface_diagnostics/checkpoint.json.
Preserve the user's poster; no remote/publication.

## Matched continuum rates and complete signed response (2026-10-05)

The full experiment-sandbox goal remains active. **446 tests pass in 234.86 s.**
Read MATCHED_MODE_SOLVER.md and MATCHED_MODE_VERIFICATION.md. New numerical
methods remain opt-in and preserve all original physical equations and strict
acceptance bounds. All new jobs are terminal. Every array and non-runtime metric
of all six signed trajectories replays bitwise. Earlier sections below are historical.

Eight separately refined continuum rate cases replay every coefficient,
physical fingerprint and non-runtime metric bitwise, with independent
quadrature 16/20. Field nr32/nz32 -> nr64/nz128 reduces radial-B/current rate
errors from 0.547444/0.480471 to 0.000244246/0.000219551. Source 7bb77514.
A scalable modal Gauss inverse and uniform physical current-inertia inner
inverse preserve all modes, physical masses and actual changing density.

Six complete 100-microsecond signed nonlinear trajectories on the explicitly
coarse nr8/nz32 / fluid4/32 mesh have three paired physical audits: 30/60 steps
at epsilon 1e-4, plus doubled amplitude at 30 steps. Source f255a576. The full
independent reference grows by 2.402175837; measured component projections at 60
steps range 2.401218278–2.402203516. Density/radial-B/current/full-E relative errors
are 0.001010704/0.001936514/0.002747908/0.008678388. Halving dt changes full E by
0.001361556; doubling epsilon changes no audited whole-domain physical response
by more than 2.07e-6. Initial representation and axis/wall/max norms are retained.
All three physical audits and two time/amplitude comparisons replay bitwise.
Stable continuous energy drift divided by the full background energy is at most
3.823e-17; that normalization does not certify small-signal energy accuracy.

Slow proposal variants and an unguarded negative-run rejection are retained.
Bounded proposals complete both signs without weakening residuals. Eleven
performance-pilot stages/endpoints have independently replayed equation audits.
A failed stale-import postprocessor has an explicit source-identity limitation;
fresh corrected audits are complete. The inspected scientific report is in
outputs/matched_mode_report. All 22 earlier verified closures and 32 preceding
archives remain unchanged. No remote or publication; preserve the user's poster.

NEXT: full growing-response mesh refinement and further time/residual controls,
then local current/enthalpy/Poynting interfaces, heated moving exact vacuum,
transport/ionization, actual coil/end boundaries, converged formation and the
experiment UI. Keep frc-7mm.7 and the full goal open. The lab geometry/gas/fill/
bias remains an exploration envelope; these are verification fixtures, not
laboratory formation or confinement predictions. NEXT_SESSION.md is current.

## Wall-space correction; complete contact transits (2026-10-04)

The full experimental sandbox goal remains active. **436 tests pass in
232.03 s.** Read STRONG_WALL_DESIGN.md and STRONG_WALL_VERIFICATION.md. The
new calculations, full replays and regression tests are terminal. Every
trajectory has independent quadrature 16/20 and complete bitwise replay. This
is bounded periodic verification, not formation. Earlier sections below are
checkpoint history; this newest section governs.

Independent reconstruction localizes the old 854533/s mode at the outer wall,
where normal flow and toroidal B violate the intended physical traces. A new
explicit Galerkin boundary space reduces that coarse rapid branch to 9203/s;
its actual full-RHS matrix agrees with the independently restricted old matrix
to 1.12125e-11. All physical-mode arrays, three matrices and matched continuum
rate cases replay bitwise. Every nonlinear stage and induction Krylov vector
satisfies the constraints directly; no accepted-state projection, damping,
artificial mass or energy/heat repair is used. The original raw Faraday defect
is separately measured, never mislabeled a passing original-incidence check.

Five new full 100-microsecond contact transits and six full fast waves are complete.
Time-difference ratios are 15.967079/15.991774 for wave B and 15.883252 for contact
density; contact energy drift decreases from 3.04416e-8 to 1.19159e-10 with time
refinement. Finest contact density/pressure errors are about 9.97995e-7 and
3.8e-6, full E error 5.00415e-4. Old frozen failures and successful pilot evidence
are retained. Source d7a422b0; diagnostic ded3bac4; continuum comparison c100c3ee.

**The growing-mode accuracy gate remains open.** An independently derived
nonuniform continuum mode gives decreasing rate errors with mesh refinement,
but the finest sampled radial-B/current rate errors remain 0.547444/0.480471.
Good unperturbed contact transit and fast-wave time convergence do not certify
that slow response. Current/enthalpy/Poynting interface budgets also remain open.

NEXT: refine field directions separately on the independent growing mode and
evolve small physical perturbations with mesh/time/amplitude controls. Then
complete local budgets, heated moving exact vacuum, transport/ionization, actual
coil/end boundaries, converged formation and experimental controls. Keep all
current numerical acceptance beads and frc-goal open. The user's prospective
geometry/gas/neutral fill/bias are exploration inputs, not lab recommendations.
Measured report: outputs/strong_wall_report/convergence.png and PDF. All 32 new
source archives are verified; outputs/strong_wall_diagnostics/checkpoint.json
records the complete checkpoint and twelve unfinished beads.

## Implicit wave methods verified; nonuniform-contact growth exposed (2026-10-04)

The full user goal remains active. **432 tests pass in 227.47 s.** All jobs from
this checkpoint are terminal. Read IMPLICIT_PLASMA_NUMERICS.md and
IMPLICIT_PLASMA_VERIFICATION.md; the inspected six-panel PNG/PDF is in
outputs/implicit_plasma_report. This is bounded filled-periodic verification,
not a completed formation simulator or a selected laboratory experiment.

Physical block preconditioning and full-equation midpoint/Gauss/Radau now have
eleven complete bitwise-replayed waves, independent quadrature16/20, original
strict induction checks and atomic time-step rejection. Physical B time ratios
are 3.975482, 15.967079/15.991775 and 7.975201/7.994196, respectively. Radau's
intentional temporal damping is measured, not called physical resistivity;
Gauss/midpoint weak-energy diagnostics retain their measured precision floor.
A captured chord-only contact failure is recovered with a Newton correction
under its unchanged original tolerance. All earlier verified physics sources
remain unchanged; initial failed candidate source archives are preserved.

**The long contact is still not validated.** Three Gauss controls and a Radau
control were stopped after large deviations from the exact solution developed.
Their terminal records preserve saved accepted checkpoints; three Gauss states
have independent physical field/energy audits. No 100-microsecond accurate
transit is claimed. The issue is not only time integration: a stationary
nonuniform-contact spatial Jacobian has a positive real mode 8.5453296e5 /s,
stable under finite-difference amplitude refinement, while uniform rest does
not have comparable growth. All three matrices replay bitwise. Paired species
pressure/current entropy feedback is necessary for the rapid mode in block
ablation controls; those incomplete controls are not physical evolution models.

An independently derived original-equation Fourier/Bessel operator verifies all
three uniform wave branches and also has slower growing nonuniform-contact
modes. Eight matrices replay bitwise; the first radial mode converges to
8763.749254 /s with axial refinement. Its radial space differs from production,
so this is not a matched-mode comparison or proof of a numerical/physical cause.
Read the explicit limitations before interpreting the spectra.

NEXT: **frc-7mm.7** is in progress. Reconstruct the production leading mode,
compare matched continuum subspaces and independent perturbation evolution,
measure mesh refinement and resolve any demonstrated coupling defect. Keep
local enthalpy/Poynting interface and weak-pressure precision gates open; then
repeat complete long contacts. Heated moving exact vacuum, transport/ionization,
coil/end boundaries, converged formation and experiment UI remain required.
The user's 1 m tube / 10–20 cm diameter / 0.1–100 mTorr neutral-fill / selectable
gas / flexible Helmholtz-like bias remains the exploration envelope. No remote
or publication. NEXT_SESSION.md is the current continuation entry point.

## Coupled moving plasma and strict residual refinement (2026-10-04)

The user explicitly made the full experiment-sandbox objective an active goal.
**425 tests pass in 219.19 seconds.** All new simulation, replay, energy,
refinement, reporting and test jobs are terminal.
The actual regular fluid/current/two-species/magnetic equations now evolve
together. The bounded filled-periodic-cylinder study has **25 cases**, all
arrays/non-time metrics replayed bitwise, with independent quadrature 16/20:
16 original parent rates, seven full fast-wave periods and two separate full
stationary spectra. Source `19f68fb08f949c4526bc1347d18652fb1eeb22d720d483a2074d6a71f7f40358`.
Read REGULAR_PLASMA_NUMERICS.md and REGULAR_PLASMA_VERIFICATION.md; raw evidence
is outputs/regular_plasma_verification and the inspected six-panel PNG/PDF is
outputs/regular_plasma_diagnostics/convergence.

At 16 elements, a full physical-mass wave period has relative density/pressure
error about 5.6764e-4, velocity 1.52793e-5, B 2.35089e-6 and J 2.59599e-6, with
explicit axis/wall/max and initial representation norms. Time-refinement ratios
are 15.999437/15.999863. Finest contact instantaneous density/electron/ion-pressure
errors are 4.04333e-10/5.65129e-10/9.10035e-9; its full E error 1.82313e-5.
The complete nonlinear contact transit has **not** been run on this assembly.

Seven additional independently reconstructed energy audits replay at both
quadratures; stable conservative-change evaluation preserves weak signals.
Finest continuous energy drift divided by wave energy is -2.97506e-8, distinct
from grid drift -2.99866e-8. Six slower-branch amplitude/mesh cases retain the
finite-amplitude nonlinear floor and an explicit weakest-pressure precision
limit; not every tiny signal converges monotonically to the continuum reference.
All 17 earlier calculation source closures are unchanged.

The initial full coordinate-Jacobian attempt is preserved with its rejection.
Captured input disproved the pressure-only hypothesis: a localized near-axis
toroidal magnetic perturbation passed the inner physical residual but failed
the full coefficient backward residual. New opt-in RefinedRegularPlasma and
RefinedKernelInduction tighten only the inner solve until **every unchanged
original outer residual** passes. Both full 360-coordinate matrices now finish
and replay bitwise; independent linear-operator errors are 2.04814e-9/2.20696e-11.
One complete reference wave is bitwise unchanged. Refinement source
`bc562813c4df25bb032e693b4e59adf0782aa9e92ee29f1233a36b34a1a2cd3f`.
frc-7mm.6 tracks full moving adoption/robustness; the original frozen class is
not silently replaced. No tolerance weakening, field fitting or energy repair.

NEXT: use the opt-in refined assembly for a complete unforced nonlinear contact
transit, with a derived stiff time method or justified explicit cost and genuine
time/mesh/positive-state controls. Close local physical enthalpy/Poynting
interface budgets. Earlier fluid accuracy limits, heated moving exact vacuum,
transport/ionization, actual coil/end boundaries, converged formation and the
experimental UI remain required. Keep frc-7mm.4/.6 and frc-goal open. This is
axisymmetric verification, not a stability/confinement or formation result.
The user's provisional 1 m tube / 10–20 cm diameter / variable gas / 0.1–100 mTorr neutral
fill / Helmholtz-like bias envelope remains recorded in RESEARCH_QUESTIONS.md.
The project has no remote; work and issue state remain local. NEXT_SESSION.md
is the current continuation entry point.

## Verified regular species component and lab envelope (2026-10-04)

The user resumed work toward an experiment-design sandbox. The next species
component is now implemented and independently verified. **416 tests pass in
210.84 s**. All 24 cases (14 original pressure/entropy rates, 8 full-transit
fixed-current evolutions, 2 frozen spectra/discontinuous controls) replay arrays
and numerical metrics bitwise; independent physical quadrature uses orders
16 and 20. Source:
`45a71e5b59cbcab6bb8556d531510d68da8389b16abcc82e1e61c0bc3237b78d`.
See REGULAR_SPECIES_NUMERICS.md and REGULAR_SPECIES_VERIFICATION.md, raw evidence
in outputs/regular_species_verification and inspected report in
outputs/regular_species_diagnostics/convergence.png (also PDF).

The handoff's point-transfer proposal has growing frozen modes; naive skew
symmetrization instead has order-unity local error. Both failures are retained.
A common physical L2 scalar projection gives a consistent skew current/entropy
operator and its paired pressure adjoint on actual n,q,l,mz,Ce,Ci. It preserves
physical species masses, local temperatures, common density, grid entropy and
thermal/current work, without global energy subtraction or added heat/diffusion.
Compatible pressure gradients are separated analytically before induction.

At 32 degree-five elements, independent electron/ion entropy-rate errors are
1.30090e-6/1.23178e-6 and pressure-electric errors 9.79306e-7/1.04547e-6.
Axis/wall/max errors converge. Work backward error is <=1.32826e-16; induced
grid thermal+current+magnetic power backward error <=2.45753e-12. Degree-three
continuous electron thermal/work defect refines from -6.95401e-9 to
-2.73774e-14 W; it is not asserted identically zero. The separate original
local enthalpy audit refines to electron/ion relative defects
2.12690e-6/2.98064e-6, while retaining actual strong div(J) and the unresolved
distributional interface budget. Audit source 4c6f6653, evidence in
outputs/regular_species_enthalpy.

One full prescribed-current shear transit has signal errors 5.78249e-2,
1.87973e-3, 3.05425e-5 at 4/8/16 elements. Time ratios are 15.99985/15.99998;
both physical species and reversed current are controlled. The finest physical
quadratic drift is -3.51823e-13, separately from grid drift -1.78565e-15.
This establishes bounded fixed-current entropy transport, not nonlinear moving
plasma stability or collisionless laboratory lifetime. All 15 earlier source
closures are unchanged. A final input-dtype guard also replays all 24 numerical
cases bitwise; its initial archive is preserved. All jobs launched here are
terminal. No remote operation or publication occurred; the repository is local.

The user supplied a provisional lab envelope: approximately 1 m long quartz
tube, 10–20 cm diameter, selectable gas and 0.1–100 mTorr neutral fill pressure;
bias field is undecided, with a Helmholtz-like arrangement preferred if suitable.
Gas temperature and ionization remain independent density inputs. These are
design variables, not validated experimental recommendations. Details are in
RESEARCH_QUESTIONS.md; existing smaller numerical fixtures are separate.

NEXT under **frc-7mm.4**: assemble the verified bulk fluid, current advection/
drift stress, species pressure/entropy and magnetic terms using the same actual
continuity; verify full independent parent equations and complete unforced
wave/contact trajectories with time/mesh controls and local physical
enthalpy/Poynting interface budgets. The issue remains open, as do frc-7mm.1/.2,
frc-7mm and the full FRC goal. Earlier fluid order/spectral accuracy, heated
moving exact vacuum, regime-appropriate transport/ionization, actual coil/end
boundaries, converged formation, experiment UI and 3-D stability/confinement
remain unproved. NEXT_SESSION.md is the current continuation entry point.

## New-thread handoff (2026-10-04)

The user explicitly stopped this thread to restart in a new one. Research here
is paused; the project and its ten unfinished beads remain open. Start the new
thread from [THREAD_HANDOFF.md](THREAD_HANDOFF.md), which lists all remaining
issues, verified evidence, original scope and the concrete next implementation.
Latest research commit: `7b3c8b2`. The interrupted continuation only read files;
no new species implementation or simulation was started. The new-thread handoff
also preserves the original user brief in ORIGINAL_PROJECT_REQUEST.md.


## Resumed checkpoint: verified current transport (2026-10-04)

The user resumed research. The saved current-transport candidate now has a
complete immutable 16-case degree/resolution/quadrature/density study, original-
platform bitwise replay, and independent physical reconstruction at orders
16 and 20. Full regression: **408 tests passed in 228.47 s**. Frozen source:
`2f682eece854d54b3b0d33026709f8d34617cc4633d8b3875a0f3e2b8e0d9540`.
All 13 previous calculation source closures were checked unchanged.

Current advection pairs the formerly missing density-energy exchange; drift
stress pairs stretching work; self drift has zero current work to roundoff.
Original electric and physical force errors decrease with refinement, including
axis and wall. At 32 degree-five fluid elements, continuous force component
errors are 2.55272e-8, 1.72774e-6, 1.00102e-5. Grid total-power backward error
is <=3.04005e-14. The independent continuous energy derivative is also measured:
coarse degree-three error is 2.35289e-5 W and refines to 3.72857e-11 W; do not
confuse grid conservation with exact continuous conservation. See
CURRENT_TRANSPORT_NUMERICS.md and CURRENT_TRANSPORT_VERIFICATION.md.

The older remote source-3f nr32 contact completed all 4116 steps and passed its
original-platform bitwise final-step and independent order-18 audit. All 14
face-radial trajectories are now complete and audited. Retrieved raw/source
hashes match; evidence is outputs/face_radial_remote_completion. frc-5b5 is
closed. No job was restarted. Remote execution was CPU NumPy/SciPy, not GPU
acceleration.

frc-7mm.3 is complete as a component audit. frc-7mm.4 is now in progress:
a four-resolution species pressure/entropy diagnostic shows why nonlinear
enthalpy gradients and entropy transport must share the discrete charge-flux
contract. Its compatible-potential comparison reduces the work mismatch to
roundoff while retaining a measured, convergent local interpolation error.
Both quadrature orders replay bitwise. Probe source c6a23e1a; complete source
and result hashes are in outputs/species_charge_probe. This is not an
implemented species integrator. See REGULAR_SPECIES_COUPLING_DESIGN.md.

NEXT derive and verify that common current/entropy flux with the actual regular
state and both species energy gradients; test original local pressure/entropy
rates before complete unforced coupled trajectories. The older fluid order/
spectral limits under frc-7mm.1 remain. Exact heated moving vacuum, physical
closures/ionization, boundary models, formation and 3-D stability/confinement
remain open. No production solver was replaced. All new study/test/audit jobs
are terminal; the project remains local with no Git remote or publication.
See NEXT_SESSION.md for the concrete continuation entry point.


## Previous user-requested pause and handoff (2026-09-30)

Paused for next week to conserve credits. See [NEXT_SESSION.md](NEXT_SESSION.md)
for the exact restart point, ordered beads and remote-study observation.
The overarching goal remains open. No new research study was launched for
this handoff; no remote job was stopped or restarted.

Beads **frc-7mm.3**, **frc-7mm.4**, **frc-7mm.5** track current-transport audit,
remaining species coupling and moving-vacuum/formation validation respectively.
The saved transport candidate passes **4 focused tests in 1.66s**; it has no
full regression or frozen 16-case audit yet. Its study CLI still requires the
records/provenance module. Keep it distinct from the **402-test / 16-case**
verified coupling baseline below. No production solver was replaced.

The updated viewer, poster and camera-orbit movie are in
`outputs/frc_showcase_checkpoint/`, using frozen verified coupling case 07.
They show actual saved vector fields in a filled periodic cylinder, not a
new formation result. Camera motion is not physical time evolution.


## Last fully verified checkpoint: physical regular-fluid/current bridge (2026-09-30)

The previous goal turn was verified progress: it preserved the failed Gauss-r
representation and verified regular-axis mixed fluid transport. This turn
connects that actual point density and physical velocity to compatible current
and magnetic work. No cell-mean reinterpretation or alternate continuity was
introduced. Force work, actual integrated axial force and angular torque use
the derived physical adjoints. Physical current inertia and its density-energy
pullback now use the same actual regular density and continuity derivative.

The full suite passes 402 tests in183.03s. Sixteen degree/resolution/quadrature/
density cases pass original-platform bitwise array/metric replay and independent
root-polynomial/FFT physical integration. Frozen source is
8c9339a6eaf0c88b88bffe108d5ca2c535a8bfbf33fa4bd11b0c0aaebba6c784.
All twelve previous source closures are unchanged. All new local study, test,
audit and reporting jobs are terminal. The six-panel report was inspected.

Degree5 smooth-force errors at2/4/8/16 fluid elements decrease to radial
9.10329e-10, azimuthal4.26004e-10 and axial8.67456e-12; first-element errors
reach roundoff/quadrature limits. Force work, integrated force and torque agree
within5.75e-16. Actual Ampere/Lorentz force is separately checked: finest global
component errors7.44700e-8,2.06750e-7,1.04214e-6; first-element axial1.38035e-6.
Manufactured magnetic-rate errors at8/16/32 field cells decrease4.62633e-4,
1.58115e-5,5.05377e-7. Separate field/fluid-grid controls expose both sources
of error. Actual continuity replaces analytic n_t in changing-density inertia.

The ideal force/electric pair agrees within9.18e-16; Hall work within7.59e-17.
The current-energy density pullback matches complex differentiation within
4.36e-16. Original induction residuals remain below1e-10. Independent physical
energy/work integration agrees closely; its largest relative discrepancy is
3.68e-13 for a nearly cancelling coupled magnetic integral. See the full audit
and REGULAR_COUPLING_NUMERICS/VERIFICATION for exact scopes and norms.

IMPORTANT: this remains an instantaneous, incomplete moving pair. Its total
energy rate contains integral mu*J²*n_t/(2e²*n²). At the finest matched fixture,
measured total rate is−0.001498461177W versus derived exchange−0.001498461105W.
The missing current-advection/drift-stress transport is explicitly retained,
not relabeled heat or canceled with a force/energy repair. No complete coupled
time evolution or formation claim follows from this component's successful audit.

The older frozen remote3f study remains13/14 complete and audited. At
2026-09-30T18:06:01.201528Z the nr32 contact is live at2050/4116 accepted steps.
Parent147751 and workers147839/147845/147852 are still in the correct
/home/m/Code/frc-tool-comparison/face-radial-20260930 directory on
m@192.168.0.202. Saved status is regular_coupling_diagnostics/remote_status.json.
Never restart from an observation timeout; no GPU acceleration is claimed.

NEXT continue frc-7mm.2: derive the shared current-advection and drift-stress
work using the actual regular continuity and cylindrical vector derivatives.
The continuum cancellation identities are now explicit in REGULAR_COUPLING_NUMERICS
and a candidate weak advection identity in REGULAR_CURRENT_TRANSPORT_DESIGN;
a global fitted cancellation is insufficient. Complete physical species
pressure/entropy and enthalpy/Poynting transport, then original parent-rate
and unforced wave/contact evolution controls before adoption. Prior unforced
fluid order and degree3 spectral limits remain recorded under frc-7mm.1.
Heated moving exact vacuum, physical closure/ionization, real ends/coils,
converged formation and3-D stability/confinement remain open. No production
solver was replaced. Full goal remains open and is now paused; state is local without a remote
Git repository or publication.

## Previous checkpoint: regular-axis mixed fluid transport (2026-09-30)

This continuation found and resolved a physical representation defect in the
previous fluid prototype. Mixed radial/axial/rotating flow revealed large
between-point velocity spikes despite conserved grid energy. At degree 5,
coarse old unforced grid/continuous energy gaps reached 0.06675, while a
physical axis-metric eigenvalue reached 576461. A naive physical force bridge
also failed the analytic f(r)=r test with nonconverging relative error 0.66543.
The old solver and all fourteen comparison trajectories remain preserved.

The new regular-coordinate state uses x=r²/2 and actual n,q=n*vr/r,
l=n*vtheta/r,mz,Ce,Ci. True angular momentum is 2*x*l. Axis velocities are
dependent regular traces with no zero-energy independent axis velocity slot.
Derived geometric work and joint entropy-density fluxes conserve grid physical
energy and angular momentum. The constant-density first-element velocity
metric is also the exact continuous physical metric. No energy reset, force
fit, artificial viscosity, radius/density floor or atmosphere was introduced.

The full suite passes 398 tests in 181.18 s. Eight original-equation rate
cases, eight complete acoustic spectra and fourteen mixed trajectories have
source/raw hashes and bitwise/independent physical audits: thirty new records,
plus fourteen preserved old comparison trajectories and three metric probes.
New source is 06dfff4762b721a154d9b8522f935427da0a13f335ac74a69549a08cdc958e7e;
old mixed source is a3d0b0c387bd4696fca958c20038d41ad7725746f4809d4052f80c36493e07b6.
Ten earlier source closures and both current closures are verified unchanged.

At degree 5 and 2/4/8/16 radial elements, manufactured continuous density
errors decrease 1.82947e-5,3.62103e-7,7.05634e-9,1.64671e-10. Fixed central
velocity RMS/c errors decrease 9.69768e-6 to 4.29509e-11. Unforced sampled
speeds remain near 0.29646*c on every mesh; final grid/continuous energy gaps
refine 5.72337e-8,5.16148e-9,3.63098e-10,1.70248e-11. Independent higher
physical quadrature agrees within 2.76914e-16 across all saved new snapshots.

Important limits remain measured, not concealed: unforced successive-mesh
velocity differences converge only roughly second order, and the degree-3
third sorted positive acoustic mode is not the third Bessel root. Higher
manufactured order does not resolve either issue. Grid conservation is not
exact continuous variable-density energy/angular conservation. The prior
part of the accumulated manufactured-work ledger is not independently
replayed; the last increment and every saved physical history are audited.
See REGULAR_COORDINATE_NUMERICS/VERIFICATION and the inspected six-panel
report in outputs/regular_euler_diagnostics. All new local jobs are terminal.

The older frozen remote3f study remains 13/14 complete and audited. At
2026-09-30T17:36:47.990170Z, its nr32 contact checkpoint is 1435/4116 accepted
steps. Parent147751 and workers147839/147845/147852 remain in the original
/home/m/Code/frc-tool-comparison/face-radial-20260930 directory on
m@192.168.0.202. Status is saved in regular_euler_diagnostics/remote_status.json.
Do not restart from an observation timeout. No GPU acceleration is claimed.

NEXT frc-7mm.2: derive and independently test the current/magnetic bridge on
the actual regular-coordinate density and velocity fields, with physical
variable-current inertia, changing-density work and true angular torque.
Start with arbitrary smooth force/torque and variable-density pullbacks; do
not silently map the new point state to the old cell/face representation.
frc-7mm.1 stays open for remaining transport limits; frc-5b5 awaits the final
remote audit. Full stress and enthalpy/Poynting transport, heated moving
exact-vacuum edge, physical closure/ionization, real boundaries, converged
formation and 3-D stability/confinement remain unproved. No coupled production
solver has been replaced. This is verified progress toward an active goal.
All repository and issue state stays local; there is no Git remote/publication.

## Previous checkpoint: cylindrical energy and angular transport (2026-09-30)

This continuation makes verified component progress; the full FRC goal remains
active. A separate cylindrical Gauss/SBP point state now jointly transports
density, physical radial momentum, **actual angular momentum**, axial momentum
and both species entropy densities. Its derived centrifugal/pressure work pairs
with actual physical kinetic/thermal energy in the declared grid quadrature.
No energy reset, force fitting, viscosity, density/radius floor or atmosphere
was introduced. Existing coupled plasma solvers were not silently replaced.

The full suite passes391 tests in165.29s; six focused tests pass again after a
reporting-only correction. Seventeen full2-D original Euler rate cases, six
independent Bessel/acoustic spectra and seventeen unforced trajectories pass
source/raw hashes and original-platform bitwise/independent physical audits.
New frozen source is14223ded4171ee70b3b790323913b1ff302d994df98a13c6328522273ebcf6ef.
All nine earlier source closures are verified unchanged.

Degree5 first-element radial acceleration errors decrease3.09383e-4,8.64907e-6,
1.24677e-6,9.05808e-8,5.87185e-9 at2/4/8/16/32 radial elements. Finest maximum
pointwise relative axis error remains2.63e-7, separately exposed. Full errors
reach9.92e-10. Degree7's fixed128-axial floor is confirmed by256-point controls,
which reach1.29e-11; its finest axis errors become nonmonotonic near1e-11.
These element counts are not equivalent to earlier face-momentum mesh counts.

Grid physical-energy rate defects are <=6.68e-15; state-rate conservation
<=8.74e-15. Each radial element's energy balances its interface flux in tests.
Full-contact density errors decrease2.37e-4/4.26e-6/6.92e-8 at16/32/64 axial
points. Fixed32-point energy drift decreases−1.50030e-11/−9.30478e-13/−5.79737e-14
under step halving. History invariants conserve within3.39e-15. Actual nonlinear
radial waves complete a full period; amplitude controls distinguish physical
nonlinear departure from their linearized Bessel target. Rotation runs span one
radial sound-crossing, not a complete rotation or plasma stability experiment.

Independent continuous fields use even polynomials in r² for n,mz,Ce,Ci,mr/r,
L/r², so the axis is regular and kinetic energy finite. This differs from the
SBP differentiation polynomial in r. Grid conservation is not exact continuous
energy/angular conservation. Across saved trajectory snapshots, higher physical
quadrature agrees within4.35e-16, while grid/continuous energy differs by up to
8.83e-11. Rotating continuous angular drift refines1.40e-9→3.71e-12→1.84e-13.
All measured differences remain in the evidence, not repaired away.

The first reporting attempt is preserved in cylindrical_entropy_before_velocity_norm:
relative velocity error was undefined at a zero-velocity turning point. Correct
RMS/sound-speed normalization changes no solver array;153 shared raw files replay
bitwise. Final evidence and visually inspected six-panel figure are in
outputs/cylindrical_entropy_verification and cylindrical_entropy_diagnostics.
See CYLINDRICAL_ENTROPY_NUMERICS/VERIFICATION for equations and limitations.

The older remote3f nr16 contact has completed2060 steps, passed its original-
platform audit and been retrieved. Together with prior local/remote cases,
13 of14 are complete/audited. Its continuous density/Te errors are1.918e-7/
8.584e-7. Only nr32 contact remains live: at2026-09-30T16:47:37Z its accepted
checkpoint is820/4116, parentPID147751 and workers147839/147845/147852 remain
in /home/m/Code/frc-tool-comparison/face-radial-20260930 on m@192.168.0.202.
See outputs/face_radial_remote_controls/remote_status.json. Do not restart a
live worker from an observation timeout; no GPU acceleration is claimed.

NEXT frc-7mm.1: verify fully mixed nonlinear r/z evolution with mesh/time
controls, then derive the variable-current/magnetic bridge on this actual
common state. No unproved nodal-to-cell/face conversion may be adopted.
Full species/Maxwell/drift stresses and enthalpy/Poynting work, heated free-edge
exact vacuum, physical transport/ionization, real ends/coils, shocks/kinetic
validity and converged formation all remain open. Audit the one remaining
remote trajectory when terminal without changing its source. Work/issue state
is local, with no Git remote or publication. This is verified progress, not a
completed goal or external blocker.

## Previous checkpoint: joint physical energy transport (2026-09-30)

The previous goal turn was verified progress: common radial momentum corrected
the near-axis force mismatch and preserved the incomplete energy contract.
This continuation derives a joint mass/momentum/two-species entropy-density
flux that conserves actual physical energy in its declared nodal integration
norm. It adds no energy reset or force-only correction. This is a new point
state, not a silent reinterpretation of the existing cell/face state.

The complete suite passes **385 tests in 164.90 s**. All 12 original axial Euler
rate cases and 13 unforced full-transit/wave trajectories are independently
audited, with bitwise state/rate/last-step/history replay. Independent 8x/16x
physical integration agrees within 1.92e-16; nodal-versus-continuous physical
energy discrepancy is separately retained, at most 6.70e-14 in saved states.
The four-panel convergence/energy figure was visually inspected.

At 128 points and order6, physical-energy-rate relative error is 6.01370e-9;
electron/ion pressure-rate errors are 4.01490e-10/1.19765e-8. Independent local
H_t matches conservative energy-flux divergence within 2.13e-14; its global
rate defect is at most 1.02e-16. Actual nodal mass, three Cartesian momenta and
both entropy-density histories conserve within 1.79e-15. This is not yet a
cylindrical angular-momentum theorem.

Finest contact/wave continuous density errors are 1.09227e-9/3.15956e-10;
their density-perturbation errors are 5.26361e-9/4.48378e-9. At fixed 32 points,
contact energy drift decreases from −2.46222e-11 to −1.52514e-12 to −9.40984e-14
under step halving, approximately fourth order. Continuous physical-energy
integration confirms the same trend. Reverse/time controls, separate species
errors and initial representation differences remain in the records.
Global entropy conservation does not give exact local isentropy: the analytic
isentropic wave acquires specific-entropy errors that refine from 1.15e-6 at
16 points to about 9.96e-12 at 128. This limitation is explicitly retained.

New frozen source is
924a586a84df69c28fe1493bf73e285efcecc02e63b6825926f5be232f94170d in
outputs/two_entropy_verification. Equations, state contracts and remaining
scope are in TWO_ENTROPY_NUMERICS/VERIFICATION. Existing source closures and
the running 3f4120fa radial prototype are unchanged. All new local computation,
test, audit and report jobs are terminal. No production plasma adoption occurred.

Seven additional remote radial-face controls completed and passed their
original-platform audits, then were retrieved to outputs/face_radial_remote_controls.
Combined with the five local cases, 12 of 14 are complete/audited. Finest nr32
acoustic n/Te errors are 1.69024e-8/1.12770e-8. The nr8 quarter-step contact
retains −4.10007e-9 physical energy drift with only −2.53e-15 finite-step
remainder, confirming that separate prototype's spatial transport defect.

Remote worker state is verified at 2026-09-30T16:00:33Z: parent PID147751 and
workers147839/147845/147852 remain in the correct study directory
/home/m/Code/frc-tool-comparison/face-radial-20260930. Saved accepted checkpoints
are nr16 contact721/2060 and nr32 contact205/4116. The other assigned cases
are complete. See outputs/face_radial_remote_controls/remote_status.json.
Never restart from an observation timeout; use existing checkpoints only after
authoritative evidence that a worker stopped. No GPU speedup is claimed.

NEXT frc-7mm.1 must extend joint transport to the actual cylindrical state:
true radial kinetic/thermal metric, regular-axis force, conservative angular
momentum, geometric torque and impermeable-wall response. Verify the same full
two-dimensional original equations and radial modes before adoption. Then
derive/verify the common variable-density current/magnetic and full species
enthalpy/Poynting work; no unproved nodal-to-cell/face bridge is allowed. Audit
the two remaining remote contacts when terminal without editing their source.
Full heated-edge/exact-vacuum, physical closure/ionization, real boundaries,
converged FRC formation and experiment requirements remain open. This is
verified component progress, not completion or an external blocker.

## Previous checkpoint: common radial physical momentum (2026-09-30)

This continuation made verified component progress; the full FRC goal remains
active. The complete regression suite passes **382 tests in 169.23 s**. Nine
original two-dimensional Euler rate cases and five completed unforced local
trajectories pass source/raw-array, bitwise replay and independent physical
energy/history audits. Seven radial acoustic spectra and eight independent
polynomial metric probes are retained. The six-panel scientific figure was
visually inspected, including corrected mesh-axis labels.

Independent interior radial faces now define the actual momentum field used by
continuity, velocity and the physical kinetic metric. The old mixed cell/flux
adjoint was inconsistent; replacing its metric alone would have been singular.
The common radial-face method reduces degree5 first-annulus acceleration errors
from the old4f sequence4.956/6.330/6.651/6.732 to
1.75954e-5/1.91939e-7/3.08493e-9/1.19443e-10 atnr8/16/32/64.
Whole-domain acceleration error reaches1.32683e-9. Actual pressure, centrifugal
and advective components are separately checked near the axis. Their individual
outer-wall constrained responses are not mislabeled as unrestricted forces.

The intermediate exact radial total-work correction is preserved with its
failed convergence test. The current prototype evolves physical centrifugal
force and exact pressure/thermal work pairing, while retaining the convective
energy discrepancy. Its scaled global energy-rate defect decreases from
8.28819e-7 to3.51141e-12 under refinement but is not exactly zero. The nr8
full-transit contact completes with n/Te errors2.90525e-5/5.15852e-5, physical
energy drift−4.10015e-9, account drift4.35e-16 and finite-step remainder−7.41e-14.
Acoustic and rotation nr8/16 also complete. No energy or temperature repair,
artificial viscosity, density floor or exterior atmosphere is introduced.

Source3f4120faf1b96ee2e94453339efe215399071fd632827aa6ae459cd9c5126737 is
frozen in outputs/face_radial_rates and outputs/face_radial_pilot. The new
equations, limits and measured results are in FACE_RADIAL_NUMERICS and
FACE_RADIAL_VERIFICATION. Existing verified source closures remain unchanged.
The prototype has not replaced the coupled moving plasma solver.

Remote state checked at2026-09-30T15:34:58Z: detached parentPID147751 and
workers147839/147845/147852 are running in
/home/m/Code/frc-tool-comparison/face-radial-20260930 on m@192.168.0.202.
Nine remaining protocol cases5,8,9,10,11,12,13,1,2 are assigned there. At the
check, nr32 rotation was complete but unaudited; accepted checkpoints existed
for nr32 acoustic90/208, nr8 quarter-step contact618/2060 and nr8 reverse
contact510/1032. See outputs/face_radial_diagnostics/remote_status.json for
the observation, not a claim that all requested cases have finished.
The remote protocol uses its CPU/SciPy backend; no GPU acceleration is claimed.

Long runs now save actual accepted state/history/energy checkpoints about every
5percent. A nonlinear interrupted/restarted wave reproduces uninterrupted arrays
bitwise. Resume only after proving the worker stopped; use --resume on its
existing directory and original platform/source. Do not restart a live job or
rerun terminal failed cases merely because an observation times out.
All local numerical/test/audit/report jobs are terminal at this checkpoint.

The previous remote4f worker59269 is terminal (SSH exit255/reset); subsequent
/proc inspection found no Python worker in its original directory. Thirteen
saved cases were audited on the original host and retrieved to
outputs/radial_work_remote_recovered:12complete,1failed with its actual rejected
step replayed. Case2, nr32contact, is explicitly interrupted/missing after a
last log at75percent and7198.1s; it is not counted complete or silently rerun.
No further hours are spent on this superseded force-inaccurate formulation.

NEXT finish/audit the remaining checkpointed controls without editing their
source. frc-5b5 remains open pending those trajectories. New frc-7mm.1 tracks
joint conservative physical convective kinetic-energy/momentum transport;
derive it without the already-disproved force-only translation correction.
Then integrate full magnetic/current/Maxwell/drift and enthalpy/Poynting work.
Heated free-edge exact vacuum, physical closures/ionization, real boundaries,
converged formation and the full FRC objective remain open. Repository and
issue state are local; no Git remote or public publication is introduced.

## Previous self-consistent fluid checkpoint (2026-09-29)

The new smooth fluid components now advance actual density, all momenta/angular
momentum and both physical log entropies. True reconstructed kinetic energy and
its changing-density work remain separate from a conservative energy account.
The latest complete regression run passes **376 tests in 151.46 s**.
No old verified current/thermal/bulk source closure was modified.

Longer tests exposed a radial acoustic instability in the first conservative
pressure stencil. Baseline source 0aa417e5 is preserved, including failed contact
states. Seven spectra with original-RHS difference and matrix-exponential/time
controls prove spatial growth. Nine uniform-translation controls also prove
that a pressure-force-only exact-energy correction would add net momentum.

Separate radial thermal/kinetic work pairing removes the growing linear modes
and strongly improves rotating force balance. Its nine 2-D parent-rate cases
and independent component audit nevertheless reveal a nonconverging relative
radial acceleration error in the first annulus (about 6.7 at nr64). Thus it is
experimental and not adopted in the coupled plasma solver. See LOG_FLUID_NUMERICS
and RADIAL_WORK_NUMERICS. Candidate source 4f0c824a662e581dc78ddfbd81a2c9cc06342fb85d19d62fda44d0fb709ab3a9.

The revised nr8 contact completes a full transit at half the original step:
density/Te relative L2 2.90525e-5/5.15852e-5, physical energy drift −4.10015e-9,
conservative account drift −1.45e-16 and separated finite-step energy defect
−7.36e-14. The quarter-step run also completes: final density/Te differences from
half step are6.29e-11/8.79e-11; remaining finite-step energy defect is−2.67e-15.
The original larger-step failure remains archived. Eight local trajectory records
and nine parent rates are independently audited. The inspected six-panel figure
and LOG_FLUID_VERIFICATION report retain the local radial-force limitation.

At that checkpoint remote candidate59269 was running at
/home/m/Code/frc-tool-comparison/log-fluid-20260929. It has since stopped;
the current recovery record above supersedes that live status. Baseline15653 is terminal
(exit1: five failed contacts, nine complete cases). All14 records are audited on
their original platform and retrieved to outputs/log_fluid_remote_baseline.
All five actual rejected steps replay with unchanged accepted arrays. Their
already unstable last states expose quadrature errors up to1.05e-3, with higher
quadrature still changing the worst by1.73e-4; this warning and the original audit
rejection are preserved, not marked converged. Both audit workers50776/40535 are
terminal. All local numerical/audit/report jobs are confirmed terminal.
Original preview11395 remains available at port8778, separate from computation.
No remote exists for Git. User-approved SSH data transfer has resumed; no public
publication occurred. Mac disk pressure was relieved by removing only66.97MB
of regenerable .venv bytecode, with about410MB free at the last check.

Four additional force-term probes isolate nr64 first-annulus pressure and
kinetic-source relative errors6.44600/6.44732, while convection is2.15e-10.
Track correction as frc-5b5; balanced rotating forces can hide these compensating
errors. The five previous verification source identities remain unchanged.

The next action at that checkpoint was to finish/audit the fixed-source
comparison and resolve the near-axis radial pressure/kinetic work adjoints before
full pressure/Maxwell/drift and current/enthalpy/Poynting coupling. The heated
free edge/exact vacuum, physical closures, formation and full FRC goal stay open.

The goal is a numerically verified and physically justified simulation in which
we can vary geometry, plasma and prescribed coil waveforms to investigate FRC
formation. Passing the existing stationary-conductor tests is not completion.

Latest comparison checkpoint: **17 WarpX GPU cases audited**, with a working
isolated RTX 4090 installation and a source-level Gkeyll fit review. Retain our
existing bench and custom axisymmetric path; an external engine is not yet a
verified shortcut around the physical plasma-vacuum edge. Full-Maxwell vacuum
waves converge at order two. The tested hybrid waves retain measured time/energy
errors; no FRC formation or GPU speedup claim is made. See ENGINE_COMPARISON.
No production backend or Stage 2 gate changed.

Latest custom checkpoint: **370 passing tests (126.72 s), 14 bulk thermal
rates and 10 changing-density trajectories independently audited**, plus six
bitwise-replayed diagnostic pressure probes and 24 linear face controls.
One conservative particle flux now advances density and both physical species
log entropies under prescribed bulk flow. The initial unstable normal-face fit
was diagnosed, preserved and corrected with centered faces, without added
viscosity or density floors. Finest full-transit density/electron-temperature
errors are 8.34e-9/6.57e-9; finest compression errors are 1.37e-6/8.27e-7.
Physical work and entropy identities, finite-step drift and initial projection
errors remain distinct. See BULK_THERMAL_NUMERICS/VERIFICATION and the inspected
six-panel report. Source a74e877fe9dd2030a6957d717d87dc1a1e96c4866fdb92ba29e5b76c157fcd2b.
All scientific workers and audits are terminal. The storage-error interruption
was recovered by recomputing only five missing cases under unchanged source.

The user explicitly authorized transfer to the GPU host after automatic review
initially rejected it. A first copy of 956,639,669 bytes of source, tests, docs,
issue export, configurations and outputs completed successfully at
m@192.168.0.202:/home/m/Code/frc-tool-comparison/recovery-bulk-a74e877-20260929/.
No Git remote was configured and no research was publicly published.
The user-requested 3D showcase is complete: an offline interactive viewer and
a visually inspected 2800×1900 current/flux-density poster from source_g1_32.
All three saved times, fixed scales, synchronized rotation and export were
checked. See SHOWCASE_3D. This initialized resistive FRC fixture is separate
from the new bulk thermal component. Only its localhost preview remains live
(session 11395, port 8778); it is not a scientific worker.

Full nonlinear adoption is still incomplete. Next derive conservative momentum
and angular-momentum evolution with the true physical kinetic metric,
velocity/vorticity, pressure/Maxwell/drift stresses and complete
current/enthalpy/Poynting work, then independent parent rates and full nonlinear
wave/contact trajectories. The diagnostic pressure force converges but is not
an exact thermal-work adjoint; it has not been adopted. The earlier smooth
solver retains fast-wave pi error 5.38315%, manufactured nr32 pi-change error
70.2663%, and slow nr64 B_t error 16.8431. Previous log-thermal, variable-current,
compatible-wave and smooth-moving source identities remain unchanged.
Heated-edge/exact-vacuum, physical transport, formation and accurate interactive
experiment requirements remain open. This is progress, not goal completion.

## Completion evidence still required

1. Conservative axisymmetric compressible dynamics: mass, momentum, pressure/
   total energy, magnetic advection and resistivity, with compatible div B.
2. Verified plasma–vacuum coupling and suitable chamber/end/outer boundaries;
   no undisclosed conducting or massive numerical atmosphere.
3. A defensible transport/ionization model for the selected regime, with local
   scale/collisionality checks; Hall/two-fluid or kinetic treatment if required.
4. Independent tests of waves, shocks, magnetic forces, diffusion/heating and
   equilibrium preservation, including convergence and energy budgets.
5. Converged formation experiments with separate reversal, topology, force
   balance, lifetime and model-validity criteria. Robustness to time step, mesh,
   boundary and transport choices; no claim of stability from axisymmetry alone.
6. Scriptable and interactive experiments with complete reproducibility records
   and diagnostics that do not overstate the supported conclusions.

## Current evidence

The initial turn was progress: Stage 0/1 were implemented with 33 passing tests
and three waveforms plus a 27-case sweep. Their topology durations are not
converged and their illustrative argon parameters raise missing-physics warnings.

The first continuation implemented a conservative dynamic r-z MHD kernel that
passes independent cylindrical fast-wave, 2-D flux-advection, exact Sod shock,
resistive-decay/heating, axis-balance, divergence and temporal-order checks.
The next continuation added limited second-order reconstruction of fluid and
staggered magnetic traces, exact regular-axis moments, total-pressure source
integration and improved one-sided wall gradients. Smooth wave/advection errors
are much smaller; axis force accuracy improves, and a particular straight-field
polynomial equilibrium is held to roundoff. First-order comparison remains
selectable and is rerun on identical grids.
The bounded-reference checkpoint had **62 passing tests**. At identical resolution the
fast-wave error is 216 times smaller and the 2-D flux-advection error 129 times
smaller than the first-order reference on the reported benchmark grids.

An independent continuum Grad–Shafranov fixture now tests genuinely curved,
closed poloidal flux with reversed Bz. It starts in equilibrium, not from coil
formation. Its artificial motion decreases with refinement, but it is not held
exactly. A lower-pressure p0=0.1 Pa stress case previously failed reconstruction
on coarse meshes despite positive cell pressure. That failure is now fixed by
local admissibility reconstruction, with consistent geometric-source integration.
Cell averages and magnetic traces are unchanged by this limiter; there are no
cell floors or energy repairs. The original failure remains reproducible with
the explicit disable switch. Limiting affects 12.10% and 2.19% of accepted
cell-stages at nr=12/24, and zero at nr=48/96 in that fixture.
The default p0=1 Pa case completes across refinement after correcting the wall
derivative. Smooth temporal refinement verifies RK2; the quarter-period wave's
limiter-switching time ratios are separately retained and less regular.
For p0=1 Pa, the finest 192×384 equilibrium run has maximum artificial speed
0.00095046 times the reference Alfvén speed after 0.1 s, with flux relative L2
drift 2.1181e-5. The initial force volume RMS converges at approximately order
1.5, so this remains a measured limitation rather than exact equilibrium balance.
Fixed-grid time refinement barely changes that drift, confirming that spatial
error dominates this particular case. Mass/energy and face divergence remain
at floating-point accuracy. Generated figures were visually inspected.

Additional equilibrium stress cases down to p0=1e-5 Pa finish conservatively,
but their coarse-grid artificial motion is not physical accuracy evidence.
An independently solved magnetized double rarefaction remains positive and
conservative with unchanged uniform B. Its bulk L1 errors converge, while the
central pressure is still 114% too high at 640 axial cells. The admissibility
limiter is inactive there: this exposes underlying flux-discretization error,
not a pressure floor. It is tracked as an accuracy requirement, not hidden by
the passing positivity checks. No full-update positivity theorem is claimed.

Measured evidence is in `docs/MHD_VERIFICATION.md`,
`docs/EQUILIBRIUM_VERIFICATION.md`, `docs/POSITIVITY_VERIFICATION.md`, and their
JSON/plot outputs. The limiter derivation is `docs/POSITIVITY.md`. These remain
bounded verification solvers, not the coupled coil-driven formation experiment.

## Dynamic vacuum coupling checkpoint

That coupling checkpoint had **74 passing tests** (18.08 s).
`VacuumCoupledMHD` advances compressible gas inside a fixed, gas-filled chamber
with reflecting slip walls/ends and magnetically transparent boundaries. A
larger magnetic-only exterior satisfies the elliptic induced-current constraint;
there are no exterior density, pressure or conductivity degrees of freedom.
An adjoint magnetic force/advection pair, implicit midpoint and gas Euler Strang
splitting account for induced self-field energy, gas energy, prescribed-source
work and positive physical Ohmic heat. This response-energy budget excludes
external-field self/mutual energy and is not a coil-circuit energy estimate.
Failure leaves accepted state, time and cumulative work intact. A source time
cap is explicit because fluid CFL cannot discover arbitrary short pulses.

Independent compact-field J cross B errors decrease from 0.2707 to 0.006437 in
volume L2 as radial resolution increases 8→64, including axis cells. Passive
advection errors decrease 0.008325→0.002155→0.0005438 at 12/24/48 cells; Gaussian
resistive-diffusion errors decrease 0.007174→0.001822→0.0004573. Fixed-grid time
halving reduces momentum error by 4.010 and 4.050. A separate infinite-space
circular-loop integral distinguishes mesh error from magnetic-box error:
at 32 cells/m, zero-flux outer extents 2/4/8 m give 5.429%, 0.6791%, 0.08840%
probe error. With independent exact boundary data the 32-cell/m error is 0.03344%.

Two prescribed coils now induce current, gas acceleration and heat in this
reference. Across the reported dynamic cases the maximum relative response-
energy residual is 1.63e-14; mass change is at most 1.42e-16. Compatible face
divergence is near roundoff. **Pulse flow is not yet converged:** final peak speed
increases 0.0661→0.1473→0.1798→0.1999 m/s at 8/16/24/32 radial cells. The magnetic
boundary also matters: at 16 cells, extent 2→4 m changes final source work by
about 19%, and 3→4 m still changes it by about 3.8%. A small solver residual is
therefore not physical boundary accuracy. Figures were inspected; full inputs,
histories, source fingerprint and versions are retained with the report.

See `docs/COUPLED_NUMERICS.md`, `docs/COUPLED_VERIFICATION.md` and
`outputs/coupled_verification/metrics.json`. Earlier reports remain historical
measurements of their recorded source, not reruns with this different solver.

## Current development step

The separate fixed-material `MaterialCurrent` component now includes physical
electron inertia, spatially prescribed positive carrier density/nonnegative
resistivity, true magnetic vacuum and explicit electron drift kinetic energy.
Conditional Coulomb transport and local Hall/inertia/drift diagnostics are
integrated into reference experiment metadata. Prescribed conductivity and
moving dynamics are unchanged. The full suite has **192 passing tests (79.10 s)**.
All 35 current-component cases, an independent time audit, 43 conditional
states and 41 prescribed-edge samples are complete. Both final figures were
visually inspected; complete source/configuration/field hashes are checked.
No numerical job remains live at this checkpoint.

CURRENT_VERIFICATION.md separates continuum error, time order, physical current
energy and numerical false reversal. The unresolved de/R=.05, nr=8 screening
case produces negative Bz down to -0.0033032 of applied bias, despite a small
energy residual; the independent exact field is positive. Finer sampled fields
remove this artifact but nr=64 still has 7.3078% magnetic-field L2 error. The
inertial decay's nr=64 current/heat errors are 0.082887%/0.10486%; initial weak
surface-current error limits convergence and is retained, not hidden.

The original Ar+ reference has di/R=4.0972 and conditional electron mean-free-
path/R=0.11088. It does not justify ignoring Hall physics. Smaller physical
edge-gradient lengths raise electron inertia, carrier-drift and kinetic
requirements. The heated sharp-edge incompatibility in SURFACE_COMPATIBILITY
remains unresolved; the new linear fixed-density component does not fix it.

The next bounded coding step is the axisymmetric Hall/toroidal reference and
independent wave-dispersion verification in EXTENDED_PLASMA_PLAN.md. Its parent
species equations, current-energy force/transport terms, separate electron/ion
thermal evolution and justified edge treatment must precede complete moving
coupling. Spatial coil force/work checks remain in MOVING_COIL_PLAN.md. Actual
gas ends, mesh/shock robustness, converged formation and the experimental UI
remain required. The original low-pressure scope is retained; a denser/easier
fixture cannot silently replace it. `Experiment(stage=2)` stays gated and the
full objective remains active.

## Infinite magnetic exterior checkpoint

That checkpoint had **83 passing tests** (20.90 s). The preceding turn
made verified progress, with no external job remaining live at this checkpoint.
`vacuum_boundary="multipole"` now represents the induced magnetic field to
infinity using regular-axis decaying Legendre modes. Boundary nodal traces,
interior magnetic energy and an independently integrated exterior energy form
one symmetric matrix. Exterior nodes and modal amplitudes satisfy stationary
vacuum rows; modes have no fluid mass, current density or resistivity. The
midpoint response-energy budget includes energy beyond the field grid.
The historical fixed zero-flux boundary remains an explicit comparison.

Independent displaced toroidal-current quadrature excites both odd and even
modes. Interior flux relative L2 errors at 8/16/32/64 fluid radial cells are
0.01257, 0.0008744, 0.0001845 and 0.00004552. Interior/exterior normal-derivative
mismatch decreases 0.01653→0.003882→0.0009413→0.0002318, including matching-surface
corners. Exterior B error at fixed 32 cells drops from 32.3% with one mode to
about 0.0364% with twelve; higher modes then reach the spatial-error plateau.
A separate analytic dipole and independent all-space volume integral check
normalization, axis fields and the sign/magnitude of exterior magnetic energy.

All fixed/open pulse comparisons were rerun from current code. At 16 fluid radial
cells and 12 multipoles, moving the matching rectangle from 2 to 4 m changes
source work by **0.00206%**. The corresponding fixed-wall change is **23.646%**,
using the 4 m result as denominator for both percentages. The prior roughly 19%
figure used the 2 m fixed-wall result as denominator. Open pulse energy residuals
are at most 2.21e-14 relative to initial response energy; mass change is at most
2.83e-16 and modal stationarity residual at most 6.32e-13. Beyond-grid magnetic
energy is positive and included. No exterior fluid or state repair is introduced.

The pulse still is not spatially converged: open-boundary final peak speeds are
0.05919, 0.13227, 0.16262 and 0.18119 m/s at 8/16/24/32 cells. At 16 cells,
halving CFL changes source work by about 0.046%, larger than the measured boundary
surface effect. More modes alone cannot fix fluid or time-discretization error.
The current single-center expansion requires a matching rectangle enclosing the
entire chamber's source sphere; efficiency for long chambers remains to evaluate.

Derivation, measured results, inspected figure, complete inputs/histories and
source fingerprint are in `docs/OPEN_BOUNDARY.md`,
`docs/OPEN_BOUNDARY_VERIFICATION.md`, and
`outputs/open_boundary_verification/metrics.json`. This verifies the open
**magnetic** exterior, not a moving free plasma surface, gas outflow, transport,
magnetic shocks or FRC formation. `frc-vacuum` and the full objective stay active.

## Magnetic shock and implicit solver checkpoint

That checkpoint had **99 passing tests** (25.17 s). That continuation made
verified progress; all measurement/report jobs completed, with none left live.

An independent strong perpendicular shock now has exact Rankine–Hugoniot jumps
and steady cylindrical branches. Finite-difference momentum and invariant checks
validate the reference independently of the numerical kernels. Both actual 2-D
solvers are tested in an interior annulus before remote boundaries dominate.
The centered coupled method creates a **false Bz reversal** at nr=48, reaching
Bz/B1=-0.027288 even with a relative response-energy residual around 1e-14.
Such magnetic structure must not be interpreted as FRC topology.

The explicit `shock_diffusion` option adds a nonnegative compressive diffusivity
that scales as h² in smooth flow and h at shocks. It is off by default and
separately records numerical magnetic heat, without changing physical Ohmic
heat or putting fluid/resistivity in vacuum. No cell floors or energy repair are
introduced. With C=1, magnetic-profile L1/B1 errors at nr=24/48/96/192 are
0.330402, 0.218542, 0.106434, 0.0581267; density errors are 0.299483, 0.187281,
0.0960794, 0.046949 kg/m³; pressure errors are 1.93626, 1.06167, 0.574753,
0.190947 Pa. No false sign reversal occurs in those C=1 windows. However,
frozen B/rho error is **not monotonic**, and nr=192 still overshoots Bz/B1 to
3.6293 versus the physical shock value 3.53234. The conservative bounded
reference preserves that special B/rho invariant to roundoff. More diffusion
is not uniformly more accurate; the coefficient study retains that tradeoff.

At nr=48, moving the inner nonanalytic extension changes pressure/B profiles by
0.3302%/0.08888% in relative L2. Doubling chamber length changes profiles below
7.1e-9; a shorter 0.01 s wall-taper control gives density change 4.86e-5.
CFL halving changes the B profile by 1.43e-5, then 3.45e-6. These are measured
controls at that grid/window/time, not general free-boundary verification.
Maximum coupled response-energy residual is 1.94e-14, mass change 2.19e-16,
and face divergence 1.25e-16 T/m. Numerical magnetic heat remains explicitly
reported (0.433 J at nr=192 for the whole chamber); physical eta and heat are
zero. Entropy-producing shock dissipation need not vanish with grid spacing.

The midpoint solve now defaults to a cached-vacuum Picard iteration at tolerance
1e-13. It solves the same implicit equations and retains atomic failure. The
original direct solve remains selectable. Direct/Picard tests cover both
magnetic boundaries and independent physical/numerical diffusion choices. Their
nr=24 shock profiles agree within 9.9e-13 relative L2; measured advance time is
6.96 times faster for Picard in that fixture. This is a CPU wall-time measurement,
not a GPU or general real-time claim.

The inspected six-panel figure, 20 completed cases, complete configuration,
initial/final diagnostics, profiles, environment and source fingerprints are in
`docs/SHOCK_VERIFICATION.md`, `docs/SHOCK_NUMERICS.md`, and
`outputs/shock_verification/metrics.json`. Old reports retain historical hashes.
Next, investigate the split method's inconsistent magnetic/fluid numerical flux
and local momentum/weak-shock accuracy, while preserving the verified vacuum
and work identities. Continue the unresolved cold rarefaction, equilibrium,
skin-layer, free-interface, transport and formation requirements below.

## Conservative fluid/CT coupling checkpoint

The current full suite has **109 passing tests** (32.92 s). This continuation is
verified progress: all 22 measurement jobs completed, including 20 successful
runs and two deliberately disabled-method failure controls. No job remains live.

`ConservativeVacuumMHD` advances mass, momentum, total fluid-cell MHD energy and
upwind constrained induction together, with shared Rusanov/UCT diffusion. It
retains the exact discrete induced vacuum and infinite magnetic exterior, with
no exterior fluid or conductivity. Vacuum-side magnetic traction and Poynting
flux include interface-current work. Positive numerical diffusivity D=a h/4 is
assigned only to plasma boundary nodes and uses the symmetric stiffness current;
its diagnosed heat is separate from physical eta J². This is declared numerical
interface stabilization, not transport. Physical and numerical electric fields
enter all adjacent total-energy fluxes, including transport along the boundary
layer, with no double-counted Joule source. An independent manufactured Poynting
identity checks that a prescribed electric curl in initially uniform B changes
magnetic energy without a spurious instantaneous thermal-energy derivative.

The original strong shock fixture reaches 0.02 s at nr=24/48/96/192. Interior
frozen B/rho defects are 2.85e-6, 9.22e-8, 4.45e-13 and 9.27e-14, with no false
field reversal. B-profile L1/B1 errors are 0.309278, 0.203389, 0.104325 and
0.0746278. Cold boundary cells require blending complete conservative/CT fluxes
toward an admissible first-order Euler proposal: 0/88, 0/176, 6/372 and 154/810
Euler stages are limited. Minimum high-order weight is 0.000332 at nr=192.
No cell pressure/density floor or energy repair is used. A failed low-order
proposal still raises atomically; no unrestricted positivity theorem is claimed.
The disabled limiter fails at 0.019111 s on nr=96; disabling interface diffusion
and limiting fails at 0.0100 s on nr=24. These failures remain in the data.

The invariant improvement is not universal profile improvement: at nr=192 the
new B error 0.07463 exceeds the historical split C=1 error 0.05813, and pressure
error is 0.36880 Pa versus 0.19095 Pa. Global limiting caused by boundary cells
reduces accuracy in the whole domain, so localization is the immediate numerical
task. Shock time halving at nr=96 changes B/density profiles by 0.02976% and
pressure by 0.07051%; doubling chamber length at nr=48 changes profiles below
2.70e-8. The original initial data and independent reference were retained.

Fluid-cell energy flux accounting is at roundoff, but independently calculated
gas-plus-full-space-induced-field response energy is not forced exact. Relative
shock response errors decrease 1.996e-4, 8.248e-5, 4.050e-5, 2.873e-5. This
explicit finite-grid error is distinct from the historical split work identity.
Mass changes are at most 2.19e-16 and compatible divergence 2.09e-16 T/m across
completed cases. All quantities remain independently diagnosed.

The same prescribed-coil pulse runs at nr=8/16/32/64. Final response error
normalized to **actual source work** decreases in magnitude 2.722%, 1.108%,
0.3406%, 0.1021%. At nr=64, physical Ohmic heat is 0.082397 J and numerical
interface heat is 0.0014725 J (about 1.8% of physical heat). Source work still
changes 0.125487→0.121817 J from nr=32→64, and peak speed increases
0.06981→0.12737→0.17880→0.20744 m/s over the four grids. The driven flow remains
unconverged. No pulse stages require flux limiting. CFL halving at nr=32 changes
work by 0.001741%, and magnetic matching extent 2→4 at nr=16 changes work by
0.002050%. Removing numerical interface diffusion at nr=32 changes work by
0.6537% and peak speed by 1.4146%, so its effects must not be ignored.

Independent smooth advection errors at nr=12/24/48 are 0.05585, 0.009702 and
0.001778; Gaussian diffusion errors are 0.007897, 0.002115 and 0.0004833, with
no flux limiting. Harmonic-polynomial tests independently verify vacuum-side
interface traces and their radial convergence. Tests cover no vacuum numerical
resistance, rest/bias preservation, cold-boundary atomic failure, local energy
conservation and independent source-budget refinement. The inspected figure,
full inputs/histories/profiles/failure records and fingerprints are in
`docs/CONSERVATIVE_COUPLING.md`, `docs/CONSERVATIVE_VERIFICATION.md` and
`outputs/conservative_verification/metrics.json`.

This remains a fixed gas-filled chamber reference. Free plasma motion and
actual end conditions, transport/ionization/Hall validity, converged formation
and a dynamic experiment UI remain unfinished. Original Stage 2 stays gated.

## Local complete-update limiter checkpoint

This continuation is verified numerical progress, not completion of the FRC
objective. The suite has **117 passing tests** (46.65 s). The default conservative
solver now chooses first-order shared face/corner fluxes only around inadmissible
Euler cells and rechecks newly troubled neighbors. Metric sources and actual
boundary energy fluxes follow that same local choice. Every flagged cell has
the independently checked low-order update, giving conditional finite
termination without a hidden global fallback. No floors or thermal repairs
were introduced; failed baseline proposals still raise atomically. The previous
global convex limiter remains an explicit, freshly reproduced control.

On the unchanged original shock, nr=192 B-profile L1 error decreases from
0.0746278 globally to 0.0593295 locally, and pressure L1 from 0.368802 to
0.334502 Pa. Frozen B/rho defect remains 8.06e-14. Only the two outermost radial
cell layers are directly affected; maximum stage support is 0.8545% of cells,
with zero direct replacements in the central independent-reference window.
Local limiting occurs more often (320 stages) but leaves distant fluxes intact.
The manufactured neighbor-chain test requires four successive detection passes;
the original nr=192 case needs at most two. Mass, fluid energy-flux accounting,
compatible divergence and vacuum stationarity remain near roundoff.

This improvement exposes a separate accuracy limitation. The nr=192 local
response-energy residual is +0.0140806 J (1.40031e-4 of initial response energy)
versus -0.00288898 J globally. The local residual is also **5.22% of the net
outward fluid energy**, a more demanding transfer scale, though that boundary
transfer is not external source work. Halving CFL raises the local residual to
0.0237089 J while central B/density profiles change only 0.001610% and pressure
0.002069%. Thus the interface evolution/energy is not time-converged simply
because the central shock is. Doubling chamber half-length at nr=96 changes
central profiles below 9.01e-9; the independent response budget remains separate.

The additional nr=384 run completes to the original 0.02 s end time: central
B/density L1 decreases to 0.0306517, pressure L1 to 0.175216 Pa, frozen-ratio
defect is 9.72e-13, and no reference-window cell is directly limited. Maximum
affected support is 0.4286% of cells, with at most two replacement passes.
However, response-energy error grows to **0.0575258 J, 30.20% of net outward
fluid energy**, and numerical interface heat rises from 0.06331 to 0.09848 J
between nr=192 and 384. The interface energy budget is therefore not converged.
Minimum density reaches 0.000388047 kg/m³ without a floor; 1575 accepted steps
take 2064 s after 101 s setup on the recorded CPU environment. The large-grid
run finished; no numerical process remains pending at this checkpoint.
All 21 case records (20 completed, one deliberately disabled-limiter failure)
and both energy audits have been source-fingerprint checked, and the final
scientific figure has been visually inspected.

The same coil pulses, time/matching controls, and smooth advection/diffusion
cases reproduce historical measured source work, response energy and analytic
errors exactly, with no local limiting. The earlier 0.1021% pulse error relative
to source work and unresolved peak-flow/skin-layer limitations remain.

An instantaneous response-rate identity now separates outward fluid power,
cell magnetic-energy rate and full-space nodal magnetic work. An independent
symmetric perturbation completes the magnetic vacuum and differences response
energy. It checks this diagnostic without changing any evolved state or budget.
The best symmetric-difference errors across the five snapshots stay below
2.62e-8 W at nr=96 and 3.65e-7 W at nr=192; probing leaves final diagnostics
exactly unchanged. At the final nr=192 state, the high-order defect rate is
-204.075 W, the first-order rate -2.83969 W, and the selected local rate
-7.23180 W. Thus this wall-state discrepancy is present in the spatial update,
well above diagnostic roundoff; its cumulative and timestep-sensitive behavior
still needs a consistent interface treatment.
The derivation is in `docs/INTERFACE_ENERGY_AUDIT.md`; measured rates, perturbation
errors, full controls, limiter maps and source fingerprints are retained in
`docs/LOCAL_LIMITER_VERIFICATION.md` and `outputs/local_limiter_verification/`.
The next interface task is to use that decomposition to improve consistent
vacuum/CT energy exchange without sacrificing local conservation or adding heat
repairs. Stronger multidimensional shocks, cold rarefaction, free surfaces,
actual ends, transport/ionization/Hall validity, formation and the dynamic
experiment UI all remain in scope.

## Compatible-energy and free-boundary-reference checkpoint

Positive staggered face-square magnetic energy and a natural variational
vacuum trace are now explicit solver options, with optional SSP RK3. Pressure,
initial energy, complete-update admissibility and diagnostics use the same
selected quadrature. Boundary work uses actual selected CT corner fields;
there are no thermal repairs or exterior fluid variables. Historical defaults
remain selectable and unchanged. Independent polynomial/force/work tests cover
unequal grid spacings, nonautonomous temporal order, third-stage atomic rejection
and mixed-boundary fallback. The pulse16 and advection12 legacy checks reproduce
historical source work, energy and analytic errors exactly.

The original shock completes nr=24/48/96/192/384. At nr=192 the response residual
falls from historical 0.0140806 J to 0.000596124 J; at nr=384 it falls from
0.0575258 J to **0.00173278 J**, or 0.728851% of net outward fluid transfer.
The fine B L1/B1 remains 0.0306581 and pressure L1 0.175312 Pa; frozen-ratio
defect is 7.07e-13. No central-reference-window cell is directly limited.
The finest run finishes 1408 steps with minimum density 0.000916625 kg/m³,
3423/4224 limited Euler stages, at most 0.48489% affected cells per stage and
three replacement passes. It takes 98.9 s setup and 2792.9 s advancement.
Its energy error still grows from nr=192, so full boundary convergence is
not established by this substantial improvement.

A per-step quadratic RK audit independently predicts the entire nr=192
response error to roundoff, with maximum stage spatial defect 2.20e-11 W and
maximum per-step identity error 1.13e-12 J. The constant-bias semidiscrete work
mismatch is resolved for this method; its finite-step error is not repaired.
At nr=24, CFL halving reduces RK3 energy error by factors 8.88 and 7.85.
The fine boundary is not time-converged: half-CFL gives 0.00100662 J with RK3
and 0.00218710 J with RK2. The full first-order spatial baseline instead gives
3.97637e-7 and 4.79246e-8 J, but degrades B L1 to 0.15656.

An explicitly experimental fixed first-order boundary strip retains the
second-order interior. Recorded coarse/fine controls test that hypothesis;
its forced order reduction is counted separately from any adaptive fallback.
It is not a production default. Source-driven skin layers, boundary accuracy
and finer shock meshes must be checked before adopting it.
The nr=64 unmodified compatible pulse residual is **0.00315436%** of applied
source work, versus historical 0.102119%. Peak flow remains mesh-sensitive.
All original pulse/smooth cases complete with inactive adaptive limiting.

The spherical FRC reference independently verifies true vacuum outside,
continuous magnetic traction with no sheet current, exact force balance and
full-space energy integrals. The existing topology detector finds one O-ring
and two axis nulls; separatrix-volume error falls from 0.219889% to 0.00374111%
over 24→192 intervals. This is initialized ideal equilibrium, not formation,
a density/transport closure, stability or confinement evidence. Its figure
has been visually inspected. See FREE_BOUNDARY_REFERENCE.md.

The full suite now has **135 passing tests (49.04 s)**. All 25 primary cases
(24 complete, one retained disabled-limiter failure), two first-order controls,
six boundary-strip probes and four RK audits are finished. Their source hashes
match the frozen engine/measurement fingerprint
`5247cc3d60224324a8f1f731b164f9fd41a16a8ea70a140af70ae1a302c022c2`.
Both scientific figures were visually inspected. No numerical jobs remain.
The nr=192 strip probe preserves B L1 0.059316 and gives response residual
3.97369e-7 J, decreasing to 4.78879e-8 J at half CFL; its 766 forced cells
and 1528 affected cells are explicitly reported separately from adaptive
fallback. Source-driven and finer-boundary checks remain open for this probe.
See COMPATIBLE_ENERGY.md and COMPATIBLE_ENERGY_VERIFICATION.md; full
source-matched records are in outputs/compatible_energy_verification.
The next numerical step is to verify the fixed boundary strip with driven
pulses and finer meshes, or improve the higher-order boundary closure without
losing conservation, CT and admissibility. Free-boundary evolution, actual
ends, transport/ionization/Hall validity, converged formation and the dynamic
experiment UI remain required. The Stage 2 formation gate stays closed.

## Boundary source, skin and end validation checkpoint

The current suite has **142 passing tests (57.97 s)**. A new independent
infinite-cylinder Bessel solution checks sinusoidal skin diffusion, regular-axis
and vacuum matching, Maxwell/Ohm identities and time-averaged surface/volume
Ohmic power. Fifteen recorded grid/time cases compare physical-only diffusion,
the production induction stencil and the fixed first-order boundary strip.
Both production choices converge in magnetic field; the surface-current
stencil has separate, slower pointwise accuracy. At nr=128 the maximum B L2
errors normalized by applied amplitude are 0.000934414 (standard) and
0.00150768 (strip), with about 7.22 cells per analytic skin depth.
This is stationary-conductor induction verification, not evolved gas dynamics.

All twelve full coil-pulse runs and four gas-end controls are complete, with
hashed lossless conserved/primitive/CT/current snapshots at 0.04, 0.06 and
0.10 s. The standard nr=16/32/64 budgets and peak speeds reproduce the previous
compatible study exactly. At nr=128, standard/strip final peak speeds are
0.221845/0.217821 m/s, but their maximum same-grid velocity L2 difference is
**9.76713%**. The 64→128 velocity changes are **7.05911% / 10.4227%**.
Thus the strip's favorable shock energy result has a material driven-flow
accuracy cost and it remains experimental; no production default is changed.
Density/pressure perturbations and surface currents also remain resolution
sensitive. A small energy residual is not an accurate formation result.

Half timestep changes pulse velocity below 0.0000710%; moving the magnetic
matching surface changes it below 0.000333%. Moving actual reflecting gas ends
from ±1 to ±2 m at a common matching extent of 3 m changes central |z|<0.5 m
velocity by 0.8174% / 0.8372%, and whole-domain source work by about 2.1%.
Spatial and material-boundary accuracy dominate those time/matching errors.
The change of gas ends is a physical-domain sensitivity, not an open-end model.

The nr=384 original shock with the strip completes to 0.02 s: B L1/B1 is
0.03065827, pressure L1 0.1753145 Pa, frozen-ratio defect 5.11e-13 and response
residual **1.01130e-7 J**, versus 0.00173278 J for the standard compatible case.
There are 985 steps and no additional adaptive limiting, but 1534 forced and
3064 neighboring affected cells. Minimum density is 0.0496321 kg/m³ without
a floor; fluid energy-flux residual is 3.28e-15 J. This preserves interior shock
accuracy and removes the observed growth of the strip's energy error through
nr=384. A separate nr=384 time study is not claimed.

All **32 new recorded cases** (12 pulses, 15 skin cases, four end controls and
the fine shock) are complete. No numerical jobs remain. Full fields, inputs,
source/environment fingerprints and the report are in
`outputs/boundary_verification` and `docs/BOUNDARY_VERIFICATION.md`.
An independent annular polynomial checks field restriction to 4.44e-16;
all saved-field hashes and the common frozen engine fingerprint match.
The final six-panel scientific figure was visually inspected after fixing
overlapping resolution-axis labels.

The next numerical implementation must address boundary/flow accuracy and
free-interface evolution using these source/skin/current checks. Keep the strip
experimental until its driven-field accuracy is adequate. Real gas ends,
transport/ionization/Hall validity, stronger shock/cold-flow and curved-equilibrium
accuracy, converged formation and a dynamic experimental UI remain in scope.
The Stage 2 formation gate and full objective remain unchanged.

## Moving material / infinite vacuum checkpoint

That checkpoint had **155 passing tests (61.39 s)**. That continuation was
verified progress toward the active objective, not completion. All 26 evolution
cases and 27 static evaluations completed. No numerical jobs remain live.
Source/configuration/environment and all lossless field hashes are checked;
both scientific figures were visually inspected, with logarithmic tick-label
collisions corrected. The full report is MOVING_PLASMA_VERIFICATION.md and the
derivation is MOVING_PLASMA.md, with records in outputs/moving_verification.

The new `MovingPlasma` reference has actual two-dimensional material motion,
adiabatic gas only inside its evolving surface, regular a=psi/r² induction
labels, analytic magnetic/pressure virtual work, an exactly curved matching
sphere and decaying magnetic modes extending to infinity. No exterior fluid,
pressure/density floors, numerical dissipation or energy repair are introduced.
This prototype is ideal-only with prescribed uniform bias; fixed-chamber
resistive references and production defaults remain unchanged.

Independent uniform-sphere/dipole energy and homologous force errors converge
approximately at second order. At nr=32 and scale=0.9, relative energy/force
errors are 0.000724094 / 0.002616710. Translating the sphere by 0.3 m exposes
spurious net axial force -0.00796587→-0.000122393 N at nr=4→32; it refines
but precludes a claim of exact magnetic axial-momentum conservation.

For exact free gas expansion at 0.1 s, velocity relative mass L2 errors are
0.0726360, 0.0268221, 0.00924936 and 0.00302204 at nr=4/8/16/32. Finest
surface-position error is 0.000138224 m and pressure-cell-average relative L2
error 0.000101608. Actual material mass is unchanged. The exact sphere radius
is 1.009950494 m; polygonal geometry and cell entropy projection are not repaired.

For the initialized spherical FRC, initial acceleration mass RMS decreases
0.713087→0.243325→0.0828970→0.0285246 m/s². Maximum artificial speed after
0.04 s decreases 0.0764241→0.0327842→0.0120209→0.00635480 m/s; finest surface
displacement is 0.000127825 m. This is measured equilibrium drift, not formation
or stability. Moving the matching sphere 2→3 m at nr=16 changes the artificial
velocity by at most 0.000331791 m/s RMS (12.2083% relative to the small reference
drift). Quadrature 5→7 changes it below 1.65e-12 m/s RMS; modes 4→8 below
2.15e-6 m/s. A tiny vacuum residual does not erase matching/mesh error.

Velocity Verlet and source-work integration show second-order time refinement
when the explicit cap controls actual steps. In the Gaussian uniform-bias ramp,
nr=8 response errors at dt=0.001/0.0005/0.00025 s are 0.000448658,
0.000112158 and 0.0000280390 J. However driven velocity changes **47.1666%**
between nr=8 and nr=16, and peak speed rises 0.166450→0.290958→0.524906 m/s
at nr=4/8/16. The driven ideal current-layer response is explicitly unconverged.
Source-work accuracy is not a reason to interpret it as physical formation.

Vertex psi is preserved near roundoff, while interpolated interior material
flux is not exact. A prescribed nonuniform deformation has relative flux errors
0.00386486→0.000929701→0.000229928→0.0000573222 at nr=4/8/16/32. Initialized
FRC interior flux drift at nr=32 is 1.85284e-6; source nr=16 is 0.000631995.
These are numerical interpolation effects, not physical resistivity or
reconnection. The regular-axis stretch approximation and low-energy/distorted
mesh modes remain explicit. No shock capture or remapping is provided yet.

The next resistive checkpoint below implements the scalar diffusion and heat
extension against the independent spherical-decay eigenmode. This ideal report
retains its historical source, including ideal source-flow limitations.

## Resistive moving material / vacuum checkpoint

The full suite has **170 passing tests (62.11 s)**. All 13 fixed-geometry decay
cases, independent stiff-mode controls and 20 actual moving cases completed.
Full configurations, histories, snapshots, source hashes and inspected figures
are in `outputs/resistive_moving_verification`; the method and report are
`docs/RESISTIVE_MOVING_NUMERICS.md` and
`docs/RESISTIVE_MOVING_VERIFICATION.md`. Frozen production/measurement SHA256 is
`452e51b4405b8a44fc6191f6716fa4e19e5798ba367bb84a0953ce3c6703503b`.
Report, analysis and plotting sources have separate recorded hashes.

`MaterialResistivity` uses a consistent plasma-only weighted mass and the same
full-space magnetic stiffness to solve scalar physical induction and recover
Jtheta. Implicit midpoint gives positive per-cell eta J² heat and a separate
source-work/energy check. `ResistiveMovingPlasma` atomically composes two
resistive half-steps with ideal material motion, updating both flux labels and
cell entropy. Fixed mass is retained and no exterior fluid, artificial
resistivity, pressure floor or energy repair is introduced. Eta=0 follows the
unchanged ideal evolution. Failure of the second resistive stage after trial
heating/motion leaves all accepted arrays, time and accounting unchanged.

An independent continuum conducting-sphere eigenmode matches a vacuum dipole
with continuous field at the surface. Maxwell/Faraday/Ohm derivatives and
all-space quadrature independently verify its current, energy and heat. At
t=0.1 s and nr=4/8/16/32/64, magnetic-field relative volume L2 errors decrease
0.0926223→0.00523351, current errors 0.0521209→0.000206685, and local-heat
errors 0.124056→0.000458581. Q1-gradient field error is about first order;
flux, current and cell heat are about second order. Axis-field and one-sided
surface-current errors are reported separately. A second radial eigenmode,
time, tolerance, matching-radius, mode-count and quadrature controls are complete.
The geometry is held fixed for this diffusion fixture; it is not force balanced.

A raw midpoint step of an independently diagonalized decaying spatial mode
produces a false amplitude -1/3 instead of exp(-4), despite a 2.22e-16 J energy
residual. Adaptive step doubling controls field, current and local heat, accepts
the unmodified two-half solution, and retains positive heat. Tolerances 1e-3
through 1e-6 reduce absolute modal-amplitude error 3.30e-4→4.13e-6 and preserve
the tested positive decay. This is local error control, not a global error,
positivity or arbitrarily small-field topology guarantee.

The actual moving cases start from the analytic FRC and then undergo resistive
relaxation, with or without a Gaussian uniform-bias ramp. Reversed topology is
initialized, not formed. With D=eta/mu0=0.1 m²/s, the driven nr=16→32 comparison
differs by 63.821% in mass-weighted velocity, 29.357% in magnetic change,
3.67623% in total current and 5.4924% in local heat density. Density change
differs by 54.0889%. These material-coordinate differences use coincident
material sampling and actual volume restriction, independently checked with
polynomial and geometric integration tests. They are not absolute solution
errors. Higher physical D=1 still gives 53.915% velocity difference; changing D
is a physical change, not numerical stabilization. The large flow error remains.

At nr=8, halving outer time caps gives approximately fourfold reduction in
velocity differences and response residuals. Driven energy residuals at
dt=0.001/0.0005/0.00025 s are 4.39132e-4/1.09776e-4/2.74435e-5 J; the last
velocity comparison differs by 0.0149697%. At nr=16, moving the magnetic
matching radius 2→3 m changes driven velocity by 0.539718%; mode/quadrature
changes are smaller. Tightening diffusion tolerance 1e-4→1e-5 leaves identical
accepted stages in these D=0.1 cases and therefore identical results. D=1 cases
do reject/subdivide diffusion trials. Mass is unchanged; maximum sampled
resistive-stage energy residual is 1.679e-14 J and vacuum relative residual
4.451e-13. These are distinct from full moving-state accuracy.

Next inspect spatial modes and surface thermal/force accuracy and add circular
coils at moving coordinates with the geometric external-source derivative.
An independent exact dipole/loop force fixture is specified in MOVING_COIL_PLAN.md.
Shock capture/remapping, material distortion, physical transport/ionization/Hall
regime, actual ends, converged formation and the experimental UI remain required.
All overarching issues and the Stage 2 entry remain open.

## Surface-resolution and thermal-compatibility checkpoint

This continuation implements optional `radial_grading` in `MovingDomain` and
records actual angular resolution in the resistive case contract. The default
uniform implementation is preserved: all archived arrays in the new source
nr=8/16/32 baselines are bitwise identical to the preceding checkpoint. The
frozen production/shared-measurement fingerprint is
`2749bf18bff1bdab93b1c338b7944ec7a6773d23c39076700f83a4faec9b90c8`.
All 21 dynamic and 13 independent cases completed, including the supplemental
finest quarter-time control (798 accepted motion steps). Analysis/render and
supplemental reproduction helpers are separately hashed. Full configurations,
lossless snapshots, histories, comparisons and figures are retained in
`outputs/surface_verification`; see SURFACE_NUMERICS and SURFACE_VERIFICATION.

Nine new independent tests bring the suite to **179 passing tests (80.30 s)**.
Revolved-polygon moments check graded geometry; exact free-gas expansion and
conducting-sphere decay check the actual kernels. Current/local heat remain
approximately second order under refinement. Grade-2 nr=32 current and local
heat relative errors are 0.000829269 and 0.00182061; Q1-gradient B error is
0.0105070. Grading trades interior resolution for surface resolution and does
not automatically improve the smooth free-gas reference.

Material-field sampling now inverts the actual initial sector geometry instead
of comparing equal normalized coordinates on different polygons. Exact convex
cell intersections resolve discontinuous thermal fields before positive Gaussian
integration. Independent affine field and geometric/indicator integrals verify
this diagnostic. Across all 18 comparisons, increasing integration order 4→6
changes reported relative norms by at most 3.38e-9 of the norm. This removes
the observed diagnostic-quadrature uncertainty without smoothing evolved fields.
The velocity norm integrates actual Q1 fields with initial physical mass;
thermal fields are compared pointwise rather than conservatively averaged onto
coarse cells. These definitions differ from the preceding report and are explicit.

The new velocity differences are 34.1273% for uniform nr=16/32 and 55.2917%
for grade-2 nr=16/32. Fine grade-2 angular doubling changes velocity 15.8362%.
Successive fine time controls change it 5.39874%, then 2.75800%; the ratio near
1.96 and nonmonotonic current differences do not establish asymptotic order.
Boundary-speed differences remain 0.04457785 then 0.01810780 m/s. The finest
quarter-time run has peak speed 0.5668651 m/s, minimum density 0.0555440 kg/m³,
physical heat 0.2738125 J and response residual 1.97602e-6 J. No floor, energy
repair or change of physical D is used. Spatial and time uncertainty remain.

Saved-state force decomposition exposes cancellation between large pressure
and magnetic surface accelerations. Recovered weak-current J cross B is a
separate discrete diagnostic, not an exact continuum reference. Mass remains
unchanged; maximum sampled resistive-stage energy and vacuum residuals are
1.386e-13 J and 7.566e-13. These do not establish local boundary accuracy.

The initial thermal-compatibility audit independently recovers the boundary
current from one-sided Maxwell derivatives with relative error below 4.7e-11.
The pressure equation gives positive equatorial initial p-dot of 3.75 Pa/s for
the unforced FRC and 1.379547904 Pa/s for the source case; a persistent smooth
p=0 trace under continuous-B/no-extra-stress assumptions would require zero.
The exact conditions and primary sources are in SURFACE_COMPATIBILITY. Keep
these initialized cases as stress tests, investigate a justified edge closure
or verified nonsmooth expansion, and do not substitute arbitrary eta tapers,
removed edge heat or easier fixtures for physical formation validation.

Both final scientific figures were visually inspected; colliding logarithmic
minor x-axis labels were corrected. The strict report serialization was fixed
to write diagnostic arrays as numeric lists, then the complete assembly was
rerun using all 18 verified comparison caches. All jobs are terminal. The goal
and relevant beads issues remain open with physical surface closure and the
independent spatial-coil force reference as concrete next work.

## Electron-current inertia and conditional transport checkpoint

This continuation adds `MaterialCurrent` with physical me and frozen prescribed
ne/eta profiles, counting current kinetic energy alongside magnetic energy and
positive per-cell collisional heating. Its sparse mixed Ampere/Ohm solve keeps
vacuum elliptic with no exterior carriers. Increment formulation preserves
stationary current; diagonal algebraic scaling improves axis-conditioned
roundoff without changing physics. A final exact-ideal-limit regression also
preserves tiny material flux under a much larger source-field change, with its
Ampere residual explicitly checked. Independent massless-limit and variable-
coefficient physical-volume energy checks pass. This is a separate linear
component, not a hidden change to the moving plasma or a complete edge model.

Thirteen new tests bring the full suite to **192 passed (79.10 s)**. Primary NRL
collision coefficients in SI agree with independent published cgs numbers.
Coulomb logs are never floored; invalid classical closure is explicitly
unavailable. Reference Stage 0/1 metadata now carries conditional thermal,
Hall/electron-inertia/Debye and drift diagnostics, retaining the original
prescribed conductivity and distinguishing inferred Drude rates. No absence of
a warning establishes validity; neutrals/ionization and tensors remain open.

All 35 current cases completed, including the added thin-layer nr=64 control.
Uniform inertial sphere decay, collisionless modified-Bessel screening and
an exact variable-density/resistivity profile independently test continuum
fields/current/heat. For de/R=.2, nr=8/16/32/64 current errors are
0.0249752/0.00768733/0.00247993/0.000828866; finest B and local-heat errors are
0.00523256 and 0.00104859. Finite inertia retains the weak derivative error in
the interpolated initial field, yielding roughly order 1.5 current convergence;
surface grading alone worsens this in the recorded control. This is explicit.

Strong screening de/R=.05 exposes false reversal on nr=8: minimum Bz/bias
-0.00330320 and 40.2668% of sampled material volume negative, while the analytic
field is positive everywhere inside. The amplitude is small relative to bias
but would invalidate a sign-only topology claim. Finer recorded samples have no
negative Bz; nr=64 current error is 0.00480722 and magnetic-field error 0.0730778.
No clipping, carrier floor, anomalous resistance or heat cutoff removes the
coarse artifact. Sample positivity is not a between-sample topology guarantee.

An independently diagonalized fixed spatial mode gives temporal halving ratios
4.002221 and 4.000555. Matching-radius, angular, quadrature, mode-count and
grading controls are retained. True-sphere energy quadrature orders 48/96 agree
within 9.172e-14 relative; maximum whole-run energy residual is 8.414e-15 of
the reported continuum energy scale and maximum equation residual 2.341e-12.
These do not erase current/field error or validate nonlinear formation.

The conditional scan records 43 assumed states and 41 prescribed-density-edge
samples. The original Ar+ case gives di/R=4.09721, de/R=0.0151831 and electron
mean-free-path/R=0.110877. Conditional perpendicular momentum and Z=1 parallel
resistivities are 8.57541e-5 and 4.37521e-5 ohm m; their limits are distinct and
neither is automatically installed. Fixed charge/temperature/ionization are
inputs, not equilibrium predictions. The prescribed vanishing-density edge
and its fixed current are a diagnostic illustration, not a thermal solution.

All source/configuration/lossless archive identities and companion documents
are checked. Frozen component measurement SHA256 is
`c66fb837ef6af44768c87d3a88e624f9ec75e318fd24a4cfd62a82752725908f`;
report, plotting, time audit and transport scan are separately hashed. Both
final figures are visually inspected. Full records are in
`outputs/current_verification`, equations in CURRENT_INERTIA and
TRANSPORT_REGIME, and results in CURRENT_VERIFICATION. All jobs are terminal.
EXTENDED_PLASMA_PLAN defines the next Hall/toroidal wave reference and the
remaining current-energy, thermal and edge integration requirements. The full
goal remains active; a component milestone is not FRC formation validation.

## Linear Hall / toroidal reference checkpoint

The full suite now has **208 passing tests (78.87 s)**, including 16 independent
Hall/finite-mass tests. `hall_parameters`, `hall_modes` and `hall_cylinder` add a
bounded uniform-density incompressible reference with regular poloidal/toroidal
fields and velocities. Exact center-of-mass/reduced-mass coefficients follow
parent species momentum; drift energy, bulk/magnetic exchange and positive heat
are paired. There is no exterior gas and no claim that its restricted periodic
reference boundary is a real plasma/vacuum or perfect-conductor boundary.

All 49 cylinder cases and 16 parallel modes (three time steps each) completed.
Actual cylinder field errors decrease from approximately 5.2% at nr=8 to 0.626%
at nr=64. Both helicities and prescribed scalar resistive damping are measured;
physical finite mass recovers ion/electron cyclotron asymptotes in the parallel
limit. A poloidal perturbation generates toroidal field and flow through Hall;
the MHD control does not. Physical volume integrals independently verify all
three energies, weak current and local midpoint heat; angular momentum is
conserved in this reference. The original argon target has not been replaced
by the separately labeled proton verification fixtures.

Parallel time-halving ratios are 3.998876–3.999721. Initial interpolation in the
slow cylindrical branch excites small fast secondary modes, giving nonmonotonic
apparent projection damping at coarse time steps. An independent full-spatial
helicity exponential and additional 2048/4096/8192-step controls resolve this:
relative time errors decrease 1.58885e-7 → 3.97212e-8 → 9.93030e-9. Remaining
spatial/projection error is explicitly separate. A one-step collisionless mode
still has 21.46% phase error despite roundoff energy conservation; no energy-only
accuracy claim is made.

Maximum whole-run energy-plus-heat residual is 6.63520e-15 relative, equation
backward error 2.38842e-17 and independent physical/modal energy disagreement
5.12417e-14 of initial energy. Cylinder source fingerprint is
`868ece2de1730a51d4aaf363652e3a19139c8289e046eda3188b57a0a45b0e3f`.
Full source/config/archive/history checks pass; records and inspected scientific
figures are in `outputs/hall_verification`. Equations and findings are in
LINEAR_HALL_NUMERICS and LINEAR_HALL_VERIFICATION. All jobs are terminal.

Next implement TOROIDAL_INTEGRATION_PLAN's conservative angular momentum and
toroidal induction in the bounded compressible core, verify geometric rotation/
pinch forces and nonlinear transport, then integrate Hall work and finite-mass
thermal/current transport from the parent species equations. The current
linear reference must not be inserted as an incomplete nonlinear closure.
Heated edge physics, spatial coils, actual vacuum/ends, shocks/remapping,
kinetic ordering, converged formation and experimental UI remain required.
The full goal remains active, and Stage 2 remains guarded.

## Nonlinear compressible toroidal checkpoint

The full suite passes **224 tests (85.76 s)**. Six new production modules add a
separate compressible state/solver with toroidal field, rotation, conservative
axial angular momentum and full magnetic/centrifugal stresses. Exact annular
r² moments enter pressure recovery and rotation/field energy, including the
axis. Both first/second-order zero-toroidal RHS recover the existing poloidal
solver, whose API and historical results remain unchanged. Failed final stages
leave accepted state, time, heat and limiter accounting untouched.

Scalar toroidal resistivity uses a positive compatible tridiagonal current mass
and matched Faraday/Poynting/heat work. It exactly recovers the independently
known 8Db diffusion of g=a+b r² at the axis and wall. Axis Ez is solved rather
than incorrectly set to zero; the conducting outer wall still has Ez=0. No
energy repair, current/density floor or hidden physical resistivity exists.
Massless MHD and fully filled periodic/conducting reference boundaries are
explicit; this is not nonlinear Hall or the free-vacuum formation interface.

All 33 declared cases completed. Exact rotation, pinch and combined equilibria
have maximum spurious radial RMS speed 3.99221e-16 m/s. A finite-amplitude
rotating entropy wave tests actual nonlinear angular-momentum/energy transport;
its nz=192 reconstructed density/rotation relative volume errors are
0.000126545/0.000125193. Torsional-wave reconstructed rotation error is
0.000462097; amplitude controls independently separate quadratic feedback from
spatial error. Stronger contrast and first-order diffusion controls remain
explicit, and no single stronger-contrast run is declared converged.

Independent conducting-Bessel decay reaches nr=128 (6689 accepted steps):
reconstructed toroidal field error 3.35910e-6, first-annulus g error magnitude
1.15482e-6 and physical heat relative error 2.07328e-5. Cell-average error is
separately recorded (1.95079e-6), not confused with field reconstruction error.
Physical quadrature orders 6/10 change reported relative norms by at most
8.21825e-17. Density time-halving ratios are 4.1723/3.9905 and angular momentum
3.7153/4.3801 against a finer same-grid reference, with limiter effects explicit.
Final mass/total-energy relative changes remain within 4.52346e-15/1.13086e-15.

All cases were regenerated under the final fixture/measurement contract, which
honors nondefault gamma and checks continuum fields between cell centers.
Frozen production/measurement SHA256 is
`24ce9b8cdb121167d78ab2ee75792bb9ace362695dd852cbc2384065f4dc69f2`.
Complete source/config/archive/history checks and scientific visual inspection
are recorded in TOROIDAL_MHD_VERIFICATION and `outputs/toroidal_verification`.
All solver, measurement, test and report jobs are terminal at this checkpoint.

Next implement NONLINEAR_HALL_PLAN: compatible nonlinear Hall induction and
Poynting work on the new compressible state, using independently derived
finite-sound-speed cylindrical modes and nonlinear transfer checks. Electron
pressure/thermal and finite-mass current transport, heated-edge compatibility,
spatial coil work, actual vacuum/chamber ends, shock/remapping robustness,
kinetic ordering, converged formation and the experiment UI remain required.
The full objective remains active and the original argon target is unchanged.

## Historical completion audit at the 224-test checkpoint

| Requirement | Authoritative evidence | Status |
|---|---|---|
| Conservative dynamic core | mhd and toroidal modules, independent tests, MHD/positivity and TOROIDAL_MHD_VERIFICATION reports | Bounded nonlinear poloidal/toroidal core implemented; positive traces without cell floors; full shock/update robustness and cold-flow accuracy remain bounded by tests |
| Physical source/plasma/vacuum coupling | Fixed-chamber, boundary, moving/resistive and SURFACE_VERIFICATION reports; SURFACE_COMPATIBILITY audit; Stage 2 guard still raises | Partial: moving resistive free surface and explicit surface refinement are implemented; conditional heated-edge incompatibility, unresolved flow/current, spatial coils, mesh robustness and actual ends remain |
| Regime-appropriate transport / Hall physics | Conditional Coulomb/local regime module, finite-current-inertia component and linear cylindrical Hall/toroidal reference; LINEAR_HALL_VERIFICATION and EXTENDED_PLASMA_PLAN | Partial: linear Hall/finite-mass waves and nonlinear compressible toroidal MHD verified; nonlinear Hall coupling, chosen formation closure, moving electron/ion thermal dynamics and physical edge/kinetic treatment remain incomplete |
| Full independent numerical verification | Prior dynamic/shock/vacuum/surface/current studies; conditional Coulomb checks; 49-case cylindrical Hall/16-mode parallel and 33-case nonlinear toroidal studies with independent species, transport, time and physical-field/energy verification | Partial; 224 passing tests; known false reversal on unresolved current layer, substantial thin-layer B error, moving-flow/cold-pressure/shock robustness and physical surface closure remain |
| Converged FRC formation | Only older stationary-conductor results, unconverged topology duration | Not demonstrated |
| Reproducible experimental interface | Stage 0/1 CLI/playback; dynamic and coupled benchmark scripts, histories and metrics | Partial; no coupled dynamic experiment UI |

Beads issues `frc-goal`, `frc-vacuum`, `frc-transport`, `frc-accuracy`, and
`frc-formation` preserve the remaining scope; `frc-mhd` records the completed
bounded reference core milestone. All numerical jobs at this checkpoint are complete. The optional request for target plasma parameters is not a blocker
to numerical verification or development.

The full objective remains active until all completion evidence is present.

## Nonlinear Hall spatial checkpoint: progress, failed axis gate retained

The full suite passes **244 tests (84.28 s)**. New independently checked
finite-sound-speed Hall modes derive all three branches from compressible
momentum, pressure and Faraday; original electron/ion and Maxwell residuals,
physical energy, conducting/axis traces and limiting dispersions are verified.
The declared reduction has massless cold electrons and adiabatic ion pressure.
A separate finite-amplitude variable-density analytic snapshot verifies local
Hall current, induction and Poynting identities before judging numerical work.

The experimental HallCurl/HallInduction/HallPoynting components cancel assembled
Hall magnetic work through a common electric quadrature and an exact adjoint
of the existing magnetic energy. Both poloidal and toroidal rates are present;
there is no carrier floor, Hall Joule heat, extra fluid force or energy repair.
They are **not enabled in ToroidalMHD**, whose prior verified path is unchanged.

All 42 instantaneous spatial cases completed. Nonlinear nr=8/16/32/64/128
physical magnetic-rate errors decrease 0.259461/0.132568/0.0666923/0.0333967/
0.0167042, and current error reaches 0.00235369. However, first-annulus gdot
error approaches 0.214581 instead of zero. A shrinking physical axis volume
hides that coefficient defect in global norms. This explicitly fails the Hall
integration gate; global conservation is not accepted as physical accuracy.
Local energy-transfer defect decreases to 0.000109583 relative but is not
installed as heating or used to repair magnetic energy.

Maximum assembled Hall work ratio is 2.72060e-17. Independent archive energy
replay backward error is 2.60690e-17 using the actual summation scale. Reference
quadrature 10/14 changes cell comparisons by at most 4.75452e-16; common
quadrature 5/7 changes the nr32 physical induction error by 1.97065e-15.
All three finite-pressure branch/sign Hall-rate errors refine to 3.52–4.00%
at nr32. These are spatial rates, not time-evolved waves. Direct quadratic
Hall amplitude powers are 1.999999987 and 2.000000001; the ill-conditioned
opposite-full-field subtraction and its small-amplitude failure remain saved.

Frozen production/measurement SHA256:
`35cea1122c521983ac6808a8d8f6d27c3cf430428561e2f0de45819598a93996`.
Every declared case/config/source/archive identity passed audit; both final
scientific figures were visually inspected. All test, measurement and report
jobs are terminal. Details: COMPRESSIBLE_HALL_REFERENCE, HALL_WORK_NUMERICS,
HALL_WORK_VERIFICATION and `outputs/hall_work_verification`.

The next action is to fix regular-axis polynomial current/electric/energy
compatibility, then local Poynting transfer and stable nonlinear Hall time
evolution. The original SSPRK2 fluid method cannot simply receive an undamped
imaginary-spectrum Hall RHS. Full electron/ion thermal and finite-mass transport,
heated-edge/vacuum/real-end closure, spatial coils, shock/remapping robustness,
kinetic ordering, converged formation and interactive experiments remain open.
No formation claim or original-argon target change is made. The full goal is
active; this is progress rather than completion or an external blocker.

## Compatible nonlinear Hall checkpoint: verified progress

The full suite passes **257 tests (88.52 s)**. The corrected electric/magnetic
spaces commute under physical Faraday in x=r²/2,z, with explicit compatible
magnetic energy carried through pressure recovery, fluid stages and diagnostics.
The old cell-energy path and failed Hall operator remain reproducible controls.
There is no fitted axis coefficient, density/current floor, energy reset or
extra Hall force. A nonlinear implicit midpoint couples Hall and positive
resistivity to ideal compressible fluid half steps; failed trials are atomic.

All **54 main plus 9 supplemental cases** completed. Axis gdot error falls
0.0203290→0.00526712→0.00132551→0.000331905→0.0000830089 for nr8→128, replacing
the historical nonzero 0.214581 limit. Uniform/strong-density, zero-component,
quadrature and radial-only controls retain convergence. Global Hall work ratio
is at most 4.41337e-18; independent physical energy replay is within 7.94556e-16
of its absolute work scale. Local energy-transfer defect reaches 9.59915e-5
relative and converges near second order, but is not zero or repaired as heat.

The actual coupled solver reproduces all three finite-pressure Hall branches.
At nr32, physical magnetic errors are 3.9960%/4.4048%/3.9631%; initial field
representation already contributes 3.9807%/4.3985%/3.9631%. This is a remaining
first-order spatial field limitation, not accurate-everywhere B inferred from
small budgets. Slow/intermediate/fast velocity errors are
1.0118%/0.2128%/0.2759%. Slow duration is 0.5 radian, the others a full period.
Slow mode retained amplitude improves 0.6961→0.9679→0.9969, so coarse damping
is substantial and explicit. A 100-fold smaller amplitude hardly changes its
nr16 error. Both propagation signs, Hall-free, nonlinear amplitude and solver
tolerance controls are archived; physical and coefficient errors stay separate.

Fixed-grid time-halving magnetic error ratios are 3.94368/3.98028 for the wave
and 3.94105/4.02325 for nonlinear evolution. The finite-amplitude 1-microsecond
flow comparison against restricted nr64 has nr32 magnetic/momentum/density
differences 0.2432%/0.3791%/1.2872% of the fine reference's **change**, rather
than its large background. This does not make nr64 an exact continuum oracle.
Independent conducting Bessel decay has nr64 field/axis/physical-heat errors
1.21612e-5/9.87096e-6/1.33088e-4, with time controls separating spatial error.

All evolved histories conserve total energy/mass to at most
6.63912e-16/5.40031e-16 relative. Independent initial/final energy replay, actual
terminal times, positivity, physical heat monotonicity, source/protocol identity
and lossless archive hashes pass audit. Angular momentum uses an absolute
normalization where its mean is zero. All three scientific figures were visually
inspected. Frozen engine/reference/measurement source is
`9c14e7c757b6d3a51bc0988faff5ff704cd58e2abb8b5ebaada646de7e4d4dd8`;
supplemental protocol and report hashes are separately recorded. All numerical,
report and test jobs are terminal. Details: COMPATIBLE_HALL_NUMERICS,
COMPATIBLE_HALL_VERIFICATION and both `outputs/compatible_hall_*` directories.

### Current completion audit and next action

| Requirement | Current status |
|---|---|
| Nonlinear compressible Hall core | Implemented and independently verified in a bounded conducting periodic reference; first-order physical B representation and limited shock/update scope remain explicit |
| Electron/ion thermal and finite-mass closure | Incomplete: cold massless electrons, single represented gas, prescribed eta; physical heating partition, electron pressure/enthalpy and nonlinear drift-energy terms still required |
| Plasma-vacuum, edge, ends and sources | Earlier moving/free-vacuum references retained; heated-edge incompatibility, spatial coils, real ends and unresolved flow/mesh errors remain |
| Converged FRC formation | Not demonstrated; initialized equilibria and waves are verification fixtures |
| Interactive physical experiments | Existing Stage0/1 interface only; guarded Stage2 and dynamic experiment UI remain required |

Next implement TWO_TEMPERATURE_PLAN's **discrete conservative charge flux and
electron pressure/enthalpy pairing**, then warm electron/ion evolution and
independent collision exchange. The weak work current is not automatically a
pointwise divergence-free cell particle flux; its dual identity must be mapped
and verified before thermal advection. The electron enthalpy drift is required
alongside pressure Ohm/Poynting terms. Finite-mass convective current, drift
energy and momentum stresses must follow from the same species equations.

Return thereafter to original heated-edge/vacuum/real-end physics, spatial coil
work, shock/remapping robustness, kinetic ordering and converged formation.
No easier or colder fixture replaces the original argon target. Beads goal,
transport, vacuum, accuracy and formation issues stay open. The full objective
is **active**; this checkpoint is progress, not completion or an external block.

## Conservative electron pressure/entropy checkpoint: verified progress

The full suite passes **266 tests (88.54 s)**. `DualCurrent` now converts
toroidal magnetic work into exact locally conservative integrated charge fluxes
on full-domain nodal dual volumes. Axis, periodic seam and physical wall trace
are explicit. An independent uniform-current control shows why omitting a
nonzero wall Btheta trace inserts a false return-current sheet and reverses
the last dual-cell current. General electric work includes wall Poynting power.

`ElectronPressure` pairs electron entropy drift by -J/(ne e) with pressure-driven
induction, retaining the enthalpy/pressure energy exchange. Stable logarithmic
spatial means and a magnetic-midpoint/thermal-logarithmic temporal state preserve
both total magnetic/electron energy and discrete physical logarithmic entropy.
The latter is distinct from a convenient linear entropy proxy. Constant-density,
isothermal and adiabatic barotropic fields provide independent zero-battery
limits. Exact gradient identities reduce axis cancellation without a fitted
multiplier, state repair, artificial heat or temperature floor.

All **46 declared records** complete: 17 spatial cases, 28 successful evolutions,
and one deliberately strict nonlinear failure. Pressure-battery gdot errors
decrease 0.0252946→0.00636035→0.00159239→0.000398242→0.0000995696 for nr8→128;
first-annulus error reaches 9.79904e-5 and thermal-rate error 4.10757e-4.
Current face-integral and nondefault gamma/strong-density controls converge.
Uniform-density induction is exactly zero; the isothermal nr64 residual is
1.31973e-9 T/(m s) absolute, with no invented relative zero-reference norm.

Actual electron pressure translates at the analytic electron velocity through
one full periodic transit, with stationary bulk and unchanged Btheta. Physical
continuous pressure errors at nz16/32/64/128 are
0.167636/0.0425302/0.0106771/0.00267227 relative to the perturbation. Nodal and
between-node errors are separate; quadrature uncertainty is at most 6.16e-17.
Opposite propagation and 70% amplitude controls are retained. Fixed-grid
pressure time-halving ratios are 4.04340/4.19891 for drift and
4.00765/4.04658 for finite-field nonlinear evolution; independent DOP853 tests
separately verify second-order temporal consistency.

Starting from Btheta=0, the isolated pressure battery generates 20.0123 microjoules
of toroidal magnetic energy after 1 microsecond at nr64, with the corresponding
electron thermal-energy loss. Compared to restricted nr64, nr32 magnetic and
pressure differences are 0.1122%/0.5289% of reference changes. These are discrete
coefficient/nodal comparisons, not an exact continuum or formation claim.
The strong-density/temperature nr64 case generates 278.219 microjoules and
retains positive electron pressure. Density-scale and solver-tolerance controls
show changes far below the measured mesh errors.

Maximum charge/work/entropy-rate backward errors are
1.38766e-16/1.80985e-16/1.74751e-16. Maximum accepted-history energy drift is
5.70516e-15 relative; physical logarithmic entropy drift is 3.35616e-12 per
particle, including accumulated nonlinear tolerance over long drift runs.
The strict nr64 tolerance 2e-13 control fails at zero accepted time with residual
9.95419e-13; budgets 80/160/320 retain that stagnation and exact untouched inputs.
No failed trial is relabeled as evolved data. The default tolerance is 1e-10.

The full protocol/configuration/source/archive/history and independent physical
pressure/energy/entropy audit passes. Both scientific figures were visually
inspected. Frozen engine/reference/measurement fingerprint:
`5cccf843cfe6ebbf6005d6b87261db7e9ae9decc1277a624b2ec90de8bbc78a7`.
Report-source and figure hashes, supplemental precision controls and lossless
fields are retained in `outputs/electron_pressure_verification`. All numerical,
test and report jobs are terminal. Details: ELECTRON_PRESSURE_NUMERICS and
ELECTRON_PRESSURE_VERIFICATION.

### Remaining integration and completion audit

This new component has **fixed density and no moving ions or Hall J cross B**.
Its constant toroidal-field/density trace boundary generally permits local
pressure Ez/Poynting flux; it is not the old Hall solver's conducting Ez=0 wall.
It is therefore deliberately not inserted as an inconsistent warm closure.
The existing compressible Hall path and all its original limitations remain.

Next implement TWO_TEMPERATURE_PLAN's conservative primal-to-dual mass overlap
and bulk flux map, proving that its divergence commutes with the existing
finite-volume continuity update. Then couple dual electron entropy to changing
density and primal total pressure/ion energy with explicit local enthalpy and
Poynting work, positive species partition and unified boundary assumptions.
Warm coupled waves, collision exchange/resistive heating partition, conduction/
ionization and full finite-mass current transport remain required. No local
energy reset may stand in for that derivation.

The physical heated edge, exact vacuum, real gas ends, spatial coil work,
shock/remapping robustness, kinetic ordering, converged FRC formation and dynamic
experimental UI are still incomplete. The original low-pressure argon target
is unchanged; Stage2 remains guarded. The full goal and its transport/vacuum/
accuracy/formation issues stay active. This checkpoint is concrete progress,
not completion or an external blocker.

## Moving ideal electron/ion coupling checkpoint: verified progress

The full suite passes **283 tests (94.53 s)**. `DualMass` maps actual primal
mass and Rusanov fluxes to the electron dual grid with a locally commuting
divergence, explicit axis/wall half support and periodic seam. Independent
affine physical flux integrals and random open-boundary tests verify this.
The Q1 scalar average is the volume adjoint of the overlap, giving an explicit
electron energy representation. No separate electron density is evolved.

`AdiabaticSpeciesMHD` now evolves moving ions and electron entropy with the same
two ideal RK stages. Changing n supplies electron compression work; positive
ion pressure is recovered from conservative total thermal energy. Equal gamma
keeps total fluid dynamics identical to the existing ideal solver, verified
bit for bit at both orders. No total-energy reset, temperature floor or hidden
retry is used. Resistivity is rejected pending physical heating partition.

All **36 declared evolutions** finish: 30 main and six stronger-gradient/common-
time controls. Independent translating contact and finite-amplitude compressing
simple waves satisfy continuity, momentum, total energy and both species
equations at two gammas. The finest continuous electron-pressure errors are
0.00561259 after one contact transit and 0.000592550 after half a compression
transit at nz256. Finest reconstructed ion-pressure compression error is
0.00150679; corresponding electron/ion temperature errors are
0.000670223/0.00343809. Node and physical-field norms remain distinct.
Pa/eV values and quadrature controls are preserved in the report.

The strong contact has a **negative physical electron entropy change**. The
reported (gamma-1) delta Se/(N kB) decreases in magnitude from 5.83338e-4 to
1.09066e-4 to 2.39562e-5 for nz64/128/256, while physical pe errors decrease
0.152234→0.0417378→0.0123756. This scheme has no proven entropy inequality or
shock heating partition; exact total energy does not remove that limitation.
Strong Mach .12 compression pe errors reach 0.000738578. Positive accepted and
reconstructed species states are retained throughout all recorded runs.

Fixed-grid time ratios are 4.03370/4.50765 contact and 4.10340/4.30592 compression.
Half-CFL changes are 0.401%/4.19% of the base node errors. A radial refinement
control's apparent 5.14219e-6 pressure difference becomes 1.61383e-15 when time
steps are matched, identifying the original CFL confound. Reference-density
scale changes pressure by 1.06722e-15. Opposite translation, gamma 1.4, first-order
and nearly cold-electron controls finish. Physical-norm quadrature uncertainty
is at most 2.41243e-10 absolute in relative-L2 units.

Maximum accepted-history mass/total-energy/linear-proxy drifts are
5.40031e-16/8.87805e-16/1.12490e-15 relative. Independent physical magnetic energy,
quarter-volume species density, electron/ion energy and logarithmic entropy
replay pass. Frozen engine/reference/measurement fingerprint:
`0e37c39afe2dc22bd122cd6d21ec5e21981cfa97b7d5d59e154bd263d88c73b4`.
Main and supplemental source/configuration/archive/history audits pass; both
scientific figures were visually inspected. All jobs are terminal. See
WARM_ADIABATIC_NUMERICS and WARM_ADIABATIC_VERIFICATION.
All 36 cases were regenerated under this final source after compacting history
storage into lossless columns. Every accepted value, endpoint, physical metric
and NPZ hash matched the earlier layout exactly; both figure hashes are
unchanged. Saved collections shrink from about 143 MB to 43 MB without losing
diagnostics. The earlier c4dca9... working fingerprint is superseded.

### Remaining integration and completion audit

This is a **common-velocity ideal thermal reference**, not a coupled warm Hall
solver. Its ideal electric wall trace does not validate the full pressure
electric field along a conducting wall, even in curl-free one-dimensional
fixtures. The isolated pressure/current and cold Hall components remain separate.
The earlier 63 Hall-case source hashes were checked against 16969a9 after the
new mass-flux interface changed source identity; they were not relabeled as
new-source records. The isolated electron fingerprint is unchanged.

Next derive and implement local conservative electron-current enthalpy and
Poynting transport on this moving state with mutually consistent magnetic
boundary traces, then verify warm Hall waves and nonlinear evolution. One
candidate worth deriving is natural fixed Btheta(R)=0 with free pressure Ez(R),
which makes boundary magnetic power zero without requiring constant wall n;
it is not the previous conducting Ez=0 boundary. This candidate is not yet
implemented or verified. No mapped local energy reset may replace the flux.

Then physical collisions/heating partition, conduction/ionization and full
finite-mass current/drift energy must follow. Original argon validity, heated
edge, exact vacuum, real gas ends, spatial coils, shock/remapping/entropy
robustness, converged FRC formation and dynamic experimental UI remain required.
Stage2 stays guarded; the full goal and transport/vacuum/accuracy/formation
issues remain active. This is progress, not completion or an external blocker.

## Warm Hall coupling and pressure-split limitation checkpoint

The full production suite passes **293 tests (94.10 s)**. A zero-Btheta/free-Ez
magnetic trace now supports full pressure gradients and varying wall density,
with exact curl-gradient/dual-current work identities. The moving warm reference
adds conservative pressure/Hall Poynting power and electron-current enthalpy
face fluxes. No local total-energy reset, pressure floor or hidden retry is used.
The five-stage coupling is experimental and retains eta=0.

All **58 main cases** finish and are audited. Conservative local pressure-energy
rate errors decrease .0940091→.000395013 over nr8→128; local ion-work residuals
also decrease. Warm branch2 electron/ion pressure errors at nr32 are
.00347082/.0173594, while physical magnetic error remains .0394869. The slow
branch is still inaccurate: finest electron/ion pressure errors are
.164298/.829787 and ion-temperature error is 2.34858 of its small perturbation.
Half-CFL and tenfold amplitude reduction do not remove the slow error.
Fixed-grid pe time-halving ratios are about four for isolated pressure, fast
warm waves and nonlinear coupled evolution. Nonlinear nr16 versus nr32 physical
pe/B differences are 4.19%/9.99% of the reference trajectory changes, not a
converged formation claim. Maximum main mass/total-energy/proxy drifts are
5.41e-16/2.15e-15/5.82e-15. Independent physical energy, density, species and
field-norm replay pass; all three scientific figures were visually inspected.

The boundary resistive mode retains a **-15.21% instantaneous first-annulus
g-rate error**, even though its evolved coefficient errors decrease. Evolved
time controls show finite time error too. This is a separate unresolved axis
consistency issue; resistivity is not enabled in the warm model.

Five supplemental runs finish, but the sixth exposes a stronger limitation.
At nr64 the maximum local pressure-only ion-energy error grows to .0921 Pa
despite a decreasing volume RMS; finer time steps retain about .0914 Pa.
At nr128, three dedicated replays with 320/640/1280 planned steps reject
negative ion pressure at about .95 microseconds, at cell [0,100]. Electron
pressure stays above 51 Pa, nonlinear residuals are below 1e-10 and global
energy drift is at most 2.08e-16. Accepted and rejected states are separately
archived; none of these failed trajectories is labeled completed.

An independent local continuum pressure-only symbol has discriminant
-7.40906e7 m²/s² at an interior point of the original nonbarotropic reference,
including background-current advection. Its frozen short-wave growth is
4303.79 |k| per second when |k| is in inverse meters. This rules out assuming
general hyperbolicity of that isolated split. It is not a proof that full warm
Hall is illposed or that this interior symbol caused the observed axis failure.
See WARM_PRESSURE_SPLIT_LIMITATION for the derivation and preserved evidence.

Frozen engine/reference/measurement fingerprint:
`ff529e0fab60fb0888806589056142c7040b3c66fcff0cfde669635e1fbb7c20`.
Complete data consist of 58 main completions, five supplemental completions and
three rejected replays, plus the analytic symbol scan. Main/supplemental source,
configuration, archive, accepted-history and failure audits pass. Reports hash
their helper sources and inspected figures. Earlier warm-ideal source is
unchanged; old electron records were source-verified at419c951 and old Hall
records at16969a9. They remain historical measurements of those sources.
All numerical/test/report jobs are terminal at this checkpoint.

### Immediate continuation and remaining objective

First analyze the combined nonideal principal behavior, then derive and verify
a joint Hall-plus-pressure stage with one magnetic midpoint and thermal mean,
shared current work, conservative face energy and atomic validation. Joint
fixed-density splitting is a candidate, not a proven cure: degeneracy for
k perpendicular to B may require retaining ion/density or finite-mass dynamics.
Preserve the exact failed fixture and time/axis/maximum-norm controls; do not
add a convenient bias, arbitrary heat or numerical atmosphere to hide failure.
Address the slow-wave and standalone resistive-axis errors explicitly.

Then physical collision/heating partition, conduction/ionization and full
finite-mass current transport remain. The original low-pressure argon target,
heated edge/true vacuum, actual gas ends, spatial coil work, shock/remapping/
entropy robustness, kinetic validity, converged formation and interactive
experiments are unfinished. Stage2 stays guarded. The full goal and its
transport/accuracy/vacuum/formation issues remain active. This checkpoint is
progress with a diagnosed limitation, not completion or an external blocker.

## Finite-mass current and precision checkpoint

The separate `CurrentThermal` subsystem now evolves magnetic/current energy
and both species thermal entropies jointly at fixed positive density and
constrained zero barycentric flow. Physical reduced electron/ion mass and full
nonlinear convective inertia are retained; no effective mass, heating repair,
atmosphere or state floor is added. CURRENT_THERMAL_NUMERICS gives parent species,
constraint reaction, boundaries and the full physical energy flux.

All 41 baseline cases are audited: 40 finish and the massless nr128 original
fixture rejects at its 0.75 microsecond trial, with positive ion pressure at
the last accepted state. Physical mass reduces the unchanged original fixture's
local pressure-defect RMS from 0.0591249 to 0.000254304 Pa across nr8→128;
maximum/axis errors also decrease. A finer-time control confirms this particular
error is spatially dominated. Full component/conservative energy drift stays
below 5.80e-16/3.87e-16. Exact nonlinear waves and full face powers refine, but
physical B/J interpolation errors remain roughly first order. A local massless
joint/moving symbol retains nonhyperbolic directions; no global theorem or
claim that it caused the numerical axis failure is made.

The independent intermediate convective-inertia forcing has a nonconvergent
first-annulus error near 0.2354. The separately verified physical inertia
inversion gives full-rate L2/axis errors 1.00341e-4/6.56197e-5 at nr128;
isolated convective L2/axis errors are 2.74914e-4/6.89978e-3. The positive
magnetic-space elimination is equivalent to the original static mixed system,
verified against fresh direct solves and a separate Fourier/Chebyshev continuum
reference. The interrupted direct archives remain labeled partial.

Four uniform-bias failures were reproduced exactly from the frozen source.
Their subtraction of nearly equal weak currents is now replaced by direct
Faraday-increment current evaluation. Small-signal Ohm and stored-state backward
errors are checked separately at unchanged tolerances. Eighteen ten-step
stationary controls and 15 nonlinear/weak-wave evolutions complete and are
independently audited. Maximum stationary B error/bias is 3.54e-14; tiny waves
retain the same physical refinement as ordinary waves. Energy/rate replay
roundoff limits under strong backgrounds are explicit, not hidden by a pressure
repair. The six repeated baseline metrics agree near roundoff.

Baseline source `6c3e50a2d595ee1e2bbe14c5a6679d8da0855ab99143fa802f3422543dec8cd8`
and exact report helpers are retained in a hashed source snapshot. Current
engine source is `80fd45bddda51871f7f36682e8d0bdc4495d0c41430191c25ef03e06ed47c4d2`.
Reports verify every relevant source/archive hash. Baseline convergence/limits
and the physical-rate plot were visually inspected. See
CURRENT_THERMAL_VERIFICATION for the detailed errors and reproduction paths.
The final full suite passes 316 tests in 103.41 s. All local numerical/audit/test
jobs are terminal; the redundant remote direct solve was explicitly terminated
and its process exit verified. Issue state is exported locally; no publication
or remote Git configuration was performed.

Next derive a genuinely moving finite-mass update with drift stress, changing
1/n current kinetic energy, bulk/current transport and velocity-gradient terms
in one conservative contract. First construct an independent manufactured
moving, varying-density reference that checks both parent species and local
energy. The fixed-n component must not be bolted into the moving split by only
changing its density matrix. Physical collision/heating partition, transport/
ionization, the original heated free edge and exact vacuum, real ends and
spatial coil work, shocks/entropy/kinetic validity, converged formation and
interactive experiments remain open. Stage2 remains guarded, the overarching
FRC goal remains active, and there is no external blocker or repository remote.

## Moving finite-mass rate and continuum energy checkpoint

The full suite passes **322 tests in 107.36 s**. `MovingCurrent` now includes
bulk B* advection, velocity curl, actual prescribed n_t/n and kinetic/mixed
scalar potentials with physical reduced mass. A positive magnetic-space solve
handles both poloidal and toroidal instantaneous rates. No density, momentum
or thermal trajectory is advanced by this component.

Independent manufactured parent species check continuity, momentum, Ohm,
Faraday, the B* identity and direct full energy flux at physical Ar+ masses,
the original n0=1e20 regime, equal masses and consistent massless limits.
Pointwise identities agree within 5.09e-15. Local energy differentiation
refines approximately fourth order; dropping mixed bulk/current transport
retains a finite 9.68227e-4 relative defect in the n0=1e18 argon fixture.

All **19 rate cases** finish and pass source/configuration/archive/physical
audits. For n0=1e18, nr8→128 B_t errors decrease 0.217882→0.0116252,
J_t 0.275820→0.00459108, E 1.54929→0.00643579 in physical volume L2.
Finest first-annulus errors are 0.000283210/0.00711149/0.00708714.
Volume B_t and axis J_t/E are approximately first order; axis B_t is second
order. Original-density and phase/mass/quadrature/tolerance controls complete.
The massless nr64 J_t volume error remains **0.403719**, with latest observed
order only 0.588; this numerical limitation is retained for investigation.

Removing density/vorticity/bulk terms separately changes the nr32 magnetic
rate by 1.89967%/0.146032%/28.0284%. Maximum full-E Faraday replay error is
1.3718e-11; the changing-density physical energy-derivative balance is at most
1.3690e-10 of component powers. These are instantaneous identities, not a
conservative complete-update energy theorem. Higher norm quadrature changes
relative L2/axis metrics by at most 1.33e-15 absolute.

The initial nr128 residual rejection and target nr32 energy-cancellation audit
failure are preserved with their exact source. Explicit linear refinement of
the original residual and stable polarization of the unchanged energy functional
resolve them without relaxing tolerances or repairing physical states. One
linear refinement is used at nr128 and in the tighter nr64 control. All final
cases share fingerprint
`d3984fb056edb1aa150d030937b9bce6c21d292b84ef80305b19b52f88af5489`;
the previous current/thermal fingerprint 80fd45b... remains unchanged.
The final source snapshot, report/helper/field hashes, inspected figure, full
test log and failed-attempt audit are retained. All jobs are terminal; no
remote Git, publication or background numerical job was introduced.

Next derive and implement a **discrete conservative moving finite-mass update**
with rho=(mi+me)n, actual continuity flux/n_t, a consistent regular velocity
map, drift-stress momentum/angular momentum, species thermal transport and full
bulk/current/Poynting/enthalpy energy flux. Verify local pressure/energy and
independent finite-mass compressible waves/manufactured moving trajectories
with spatial/time convergence and atomic failure. Do not bolt the rate onto
the old split. Physical collisions/transport/ionization, heated free edge and
exact vacuum, real ends/spatial coils, shocks/entropy/kinetic validity, converged
formation and interactive experiments remain open. Stage2 stays guarded and
the overarching goal remains active. See MOVING_CURRENT_NUMERICS/VERIFICATION.

## Conservative moving finite-mass integration checkpoint

The full suite passes **328 tests in 120.97 s**. FiniteMassMHD now evolves
rho=(mi+me)n, all momentum components, conservative total energy, psi/g and
electron entropy through four complete unsplit RK4 derivatives. The actual
mass flux supplies dual n_t, induction retains physical reduced mass and bulk
velocity gradients, and shared fluxes carry drift stress, species enthalpy,
mixed current energy and full computed Poynting power. All stages validate
positive species states before atomic commit. No physical step retry, thermal
floor, artificial mass, atmosphere or total-energy reset is introduced.

All **21 complete-rate cases and 18 declared trajectories** finish and pass
source/configuration/archive/physical audits. Parent manufactured species and
independent exact-mass compressible waves verify rates and final states. The
manufactured and unforced nonlinear runs last 50 ns; the fast waves last one
quarter period (29.1294 ns). Slow/medium branches have rate evidence only.
Global mass/energy-budget errors are at most 4.44e-16/3.26e-16; independent
physical energy replay differs by at most 2.55e-14. Explicit source work matches
independent parent total-energy change within 5.08e-12 J. All accepted species
pressures are positive. These conservation properties do not prove local accuracy.

The finest manufactured nr128 mass/total-energy rate errors are
.000861842/.000529452, but pe/pi rate errors remain .00774879/.143364.
The pi-rate maximum increases slightly from 30535.9 to 31209.6 Pa/s at nr64→128,
despite decreasing volume/axis errors. The slow-mode nr64 pi-rate error is
41.4570 times its physical signal. Tenfold amplitude reduction leaves these
errors essentially unchanged. Preserve maximum norms as well as L2 evidence.

At nr32, evolved manufactured ion pressure has .00362760 relative error against
the full pressure and **.783580 against the true pressure change**. The original
n0=1e20 case retains .00250372/.540816. Fast-wave magnetic/electron/ion pressure
perturbation errors are .0398052/.0105524/.595516. The weak-wave control preserves
that limitation. Fixed-grid dt halving has a much smaller effect than mesh error.
Magnetic time ratios are about 15.4; pressure ratios range 3.15–5.19, so formal
RK4 order is not asserted for every quantity. Unforced electron physical log
entropy decreases by 5.13e-6→3.95e-7 per particle under refinement; the method
has no entropy inequality or physical shock partition.

The slow-wave continuity error is chiefly the limited central reconstruction,
not Rusanov density diffusion. At nr32 those error/contribution norms are
28.6816 and 1.03981 times the true rate. A separately tested cubic reconstruction
of exact annular linear momentum averages reduces the density-rate error to
.0202495 at nr32 and .000179355 at nr128. With actual finite-amplitude momentum
it instead levels near .00311608 against the linear oracle. This candidate is
**diagnostic only**, not adopted into the coupled solver. Correcting mass flux
alone would leave enthalpy/kinetic/momentum work inconsistent.

Baseline rate source is 8937eefb61c197962881e7d844bee6d9f5dfce45aff93025c1c35e00dba63ea7;
expanded trajectory source is 3edae52928f5cdaa20ca240f11d6fa729a137be327c9b1ffc64f09c6440f42fb.
Baseline production sources are unchanged; only provenance coverage expands.
Prior moving-rate d3984fb... and current/thermal 80fd45b... fingerprints remain
unchanged. Full records, fields, histories, physical/source audits, diagnostic
arrays, source snapshots, inspected accuracy figure and test log are retained
in outputs/finite_mass_initial_attempt and outputs/finite_mass_evolution.
See FINITE_MASS_NUMERICS/VERIFICATION for scope and complete measured errors.

Immediate next work is to derive shared higher-order particle, momentum,
enthalpy and kinetic-energy transport. Verify a translating two-temperature
contact and all finite-mass wave pressure rates, including near-solenoidal
cancellation, before adopting a reconstruction. Retain independent parent
pressure/energy, axis/max norms, exact uniform flow and atomic failure checks.
Then revisit physical collision/heating partition, conduction/ionization,
original heated free edge and exact vacuum, actual gas ends/spatial coils,
shocks/entropy/kinetic validity, converged formation and the experimental UI.
frc-7mm and the overarching goal remain open; this is progress with a diagnosed
accuracy limitation, not an external blocker or a completed formation simulator.

## Shared smooth transport and magnetic cancellation checkpoint

The full suite passes **336 tests in 118.11 s**. The separate SmoothFiniteMassMHD
reference uses shared cubic conservative face transport, mean-preserving
transverse quadrature and a consistent ratio of reconstructed dual U and n for
entropy flux. The historical limited solver is unchanged. This central smooth
method introduces no physical/artificial diffusivity, floors, hidden fallback
or atmosphere. It has no shock/entropy/positivity theorem; invalid stages reject
atomically. Exact vacuum is still an unfinished coupled-model requirement.

All **31 rates and 27 trajectories** pass the independent source/field/history,
parent-reference, physical component-energy/entropy and conservation audit.
Energy-budget residual is at most 3.36e-16; normalized integral drift is at most
4.97e-15. All measurement and audit processes are confirmed terminal. Source
dd060015dc510d312f22a53588d2809ad806de02de1dadb4db956e8887ab4ead and its snapshot
are in outputs/smooth_finite_mass_verification. Earlier finite-mass 3edae529...,
moving-current d3984fb... and current/thermal 80fd45b... identities are unchanged.
The first smooth attempt, original electric-kernel failure, serialization error
and orphan fields remain preserved. All 49 old complete records and both orphan
arrays reproduce exactly; new contact records supply the missing histories.

At nr128, manufactured pi-rate error improves from .143364 to .0268603; its
maximum now decreases to 3225.44 Pa/s. Slow nr32/64 pi rates improve from
117.248/41.4570 to .0691434/.00883534. Full fast-wave nr32 pi error improves
from .595516 to .0538315; magnetic error remains .0398052. Manufactured nr32
pi error in its small physical change remains .702663. Exact global energy
and accurate thermal rates do not establish accurate magnetic dynamics.

Three unmagnetized two-temperature contacts finish a full 100-microsecond
transit. Finest nz64 nodal pe/cell-mean pi errors are 6.63e-6/.000316986;
continuous pe/pi errors are .0000750472/.0110717 and temperature errors
.000498286/.0127339. Initial projection/representation errors, reverse flow,
limited entropy and half-dt controls are explicit. Small false B/J drift is
not monotonically convergent. Three magnetized contacts last only 50 ns and
must not be described as full-transit magnetic verification. Physical entropy
defects remain and reverse sign with flow in the contact control.

The exact separable electric-kernel representation removes a contact roundoff
rejection while preserving strict small-signal residuals and original full
equation backward errors. The full recombined Faraday error relative to its
tiny rate remains explicitly .016679 in the preserved example; there is no
claim of well-conditioned full-array reconstruction or arbitrary-kernel support.

Slow magnetic-rate errors remain 1004.98/263.069/66.8838/16.8431 at nr8→64;
slow conserved-energy-rate errors remain 106.207/26.5782 at nr32/64. Exact
velocity/current forcing controls identify incompatible lower-order field
representations, amplified by roughly 19200-fold ideal/Hall cancellation.
Neither mesh refinement alone nor fitted velocity from numerical current is an
adequate resolution. Next derive shared higher-order velocity and magnetic/
electric bases with physical Hodge/inertia quadrature and commuting curl.
First verify independent polynomial fields, Ampere/current convergence, all
wave branches and full original Ohm/Faraday/energy contracts, then repeat the
complete smooth trajectories without changing their physical inputs.

The numerical/physical requirements in frc-7mm, frc-accuracy, frc-transport,
frc-vacuum and the overarching goal remain open. Continue collision/heating
partition, conduction/ionization, heated free-edge exact vacuum, real ends/
spatial coils, shocks and kinetic validity, converged formation and the UI.
No FRC stability or particle-confinement claim follows from this checkpoint.

## Higher-order compatible field and harmonic precision checkpoint

The full suite passes **343 tests in 118.34 s**. New `compatible_polynomials`
and `compatible_geometry` supply local degree-1/3/5 nodal/interval-average bases
with commuting derivative, unchanged Faraday incidence and physical tensor
electric/magnetic masses. `compatible_velocity` reconstructs actual annular
conservative averages. `compatible_current` uses physical finite species masses
and constant-background inertia. Five independent polynomial/integral/current/
energy/Faraday/mixed-system tests verify the component, not just its algebra.

The first 25-case rounded-total study exposed a precision limit: slow nr32/64
B_t errors .0117924/.0245827 despite far smaller spatial error in the weak
perturbation-only control. A hundredfold weaker full-bias wave gave 1.17396.
This complete attempt is preserved under source 485edaa1fdb9bffcb76b2a4e9196e3e4b1ddc071b20f319b4c3062b42ebae667.
An independent three-way probe distinguishes current-operator cancellation from
digits already lost in storage; subtracting bias afterwards does not recover them.

`HarmonicMagnetic` now stores an explicit constant physical bias and independently
retained increment. Uniform axial bias has exactly zero current under the stated
fixed-flux periodic boundary. Full field/local energy keep the bias and local
cross term; its cross integral is zero by exact boundary values. Two polynomial
tests verify weak currents through amplitude 1e-12 beside that unchanged bias,
small-increment energy and independent full local energy. No signal threshold,
fitted velocity/current/mass, energy repair or numerical atmosphere is introduced.

All **25 affine cases and 25 rounded comparisons** pass the independent full
source/input/field/physical-norm/original-equation audit. Main full-bias degree5
slow B_t errors at nr8/16/32/64 are .789352/.0141969/.000190057/2.36657e-6.
Finest intermediate/fast errors are 6.68e-9/2.51e-9. Slow current-rate error
decreases to 7.73307e-5, less favorable than magnetic accuracy; axis/max errors
are retained. Hundredfold weaker signal, opposite propagation, original n1e20
and quadrature controls retain the improvement. The exactly matched earlier
linear diagnosis used the same physical inputs, without analytic forcing fields.

Largest independent physical energy discrepancy is 2.75e-14, including error
normalized by the small increment energy. Original backward errors are below
1.83e-12, but full-E Faraday replay relative to its tiny rate reaches 5.88e-10
after restoring pressure potential. Rounded export still loses current accuracy;
its worst current-relative Ampere replay error is .000129808. These distinct
measures are not mislabeled as small-signal solver residuals. All jobs are
confirmed terminal, and the scientific figure was visually inspected.

Final source is 8cfb1a2f89f7ee58fadae2e47ae9e4fbb5c49037face19da117518188ed8d720
in outputs/compatible_wave_verification; report, snapshots, failures, all arrays,
test/worker logs and checkpoint hashes are retained. The smooth nonlinear source
dd060015... and prior finite-mass 3edae529... remain unchanged. The new component
is **not yet adopted in nonlinear moving evolution**, and no trajectory energy,
entropy, shock, formation or confinement claim follows from these rates.

NEXT derive positive variable-density physical current inertia and local energy
derivatives on this field basis, with actual conservative mass-flux n_t. Derive
the common species transport/pressure-gradient adjoint and full velocity/vorticity,
stress/Poynting/current-enthalpy work before replacing the old state partition.
Retain local bias cross energy and exact physical masses. Verify full independent
parent manufactured equations and coupled wave/contact trajectories, including
small-pressure-change, axis/max, time, quadrature and atomic failure evidence.
The original heated free edge/exact vacuum, transport/ionization, real ends/coils,
shock/entropy and kinetic validity, converged formation and UI remain open.


## Variable-density physical inertia and the axis thermal limitation

This continuation made verified component progress. The final suite passes
**351 tests in123.17 s**, including eight new focused checks. All14 prescribed-
continuity manufactured cases completed and passed original-equation, source,
raw-array and independent physical local-energy audits. Every worker is terminal.

`CompatibleDensity` uses actual conservative cell-average number density and the
same linear reconstruction for supplied continuity n_t. Its physical adjoint
pairs the changing-current-energy density work with the actual cell state.
`VariableCompatibleCurrent` supplies positive variable-density inertia through
tensor products, with physical masses and the unchanged original induction
solve. Local J dot J_t, −J² n_t/(2n), and full B dot B_t (including bias cross
work) remain separate. Negative/nonfinite reconstructions or unrepresentable
inertia reject, without floors or energy repair.

The14-case degree3/5, resolution, quadrature, physical-density and flow-sign
protocol is instantaneous and explicitly manufactured. Finest degree5 B_t/J_t
relative errors are4.96985e-10/4.65041e-10. Local current-energy-rate and its
density-contribution errors are8.93213e-10/1.17057e-9. Independent root-product
field and primitive-polynomial density quadrature agree with production local
energies/rates within4.57e-15. Original backward error is below9.94e-13; checked
solver residual below1.33e-11 at1e-10 tolerance. Full-array Faraday rate-relative
error reaches5.77e-10 and remains separately exposed.

`CompatibleCharge` replaces the old diagonal-g stream with (H_g g)/(dx dz),
restoring graph divergence and pressure/current work pairing. The explicit
`CompatibleSpeciesPressure` bridge retains the old lumped nodal thermal state.
Its species work and fixed-n logarithmic entropy cancel near roundoff, but its
physical pressure E remains roughly second order and entropy rate is inaccurate.
Finest electron/ion whole-domain nodal entropy-rate errors are6.35129%/6.28264%.
The first four radial rows contain99.8072% of the electron squared error; their
local relative error grows toward39.09% with refinement rather than converging.
The interior relative error is0.09767% at nr64. This local axis defect is the
next thermal target; global conservation or shrinking error volume cannot
justify adoption. The scientific figure and radial profiles were inspected.

Final source fc106dc19b9bc208aaf7d9b36a43dc69ea49f4d49d126b4f4dfd2ecd85c59e97
is saved in outputs/variable_current_verification. The complete initial source
310809e2 and all14 cases/audit/test logs remain in
outputs/variable_current_before_range_guard. Final finite-range guards reproduce
every saved numerical array bit for bit across all14 cases. The old smooth
source dd060015 and linear compatible source8cfb1a2f remain unchanged.

NEXT derive a common physical higher-order thermal density/energy/transport
metric and its pressure-work adjoint, explicitly correcting the nonconverging
near-axis lumped entropy rate. Then derive actual conservative velocity/vorticity,
stress and complete bulk/current/enthalpy/Poynting energy transport, and repeat
independent parent-species and full wave/contact trajectories. This component
has not replaced SmoothFiniteMassMHD. Heated free edge/exact vacuum, physical
transport/ionization, real end/coil coupling, shocks/entropy, kinetic validity,
converged formation, three-dimensional stability/confinement and UI remain open.
No remote exists; work and issue state are saved locally without publication.

## Physical thermal metric and full-transit current transport checkpoint

This continuation made verified component progress. The complete suite passes
**361 tests in 121.42 s**. All 18 rate cases and nine prescribed-current
transport cases completed, with original-equation, source/raw-array and independent
physical audits. Every worker, audit and report process is terminal. The
six-panel scientific figure was visually inspected.

`EntropyScalar` supplies the physical number-weighted scalar metric on the same
reconstructed conservative density as variable-current energy. The new nodal
specific log entropy in `CompatibleLogThermal` gives positive representable
pressure, a pressure field paired with physical thermal work and explicit
density/entropy energy derivatives. Current transport conserves physical log
entropy; its exponential proxy is a separately measured approximate invariant.
No pressure floor, current/mass fitting or energy repair is introduced.

Degree-5 electron first-four-row log-rate errors are 1.02400e-3, 2.96189e-5,
5.42997e-7 and 8.84960e-9 at nr=8/16/32/64. The matched old bridge approaches
0.390908 instead of converging. Finest continuous physical proxy-rate error is
1.30658e-8 versus the old representation's 0.117630 in that same norm. Finest
electron/ion pressure-E errors are 4.37e-8/4.21e-8. Independent local energy
agrees within 6.83e-15; pressure work and declared bulk-energy exchange defects
are at most 1.85e-14/1.37e-14. These conservation measures are not accuracy claims.

Constant entropy has exactly zero stored current-driven entropy/magnetic rate.
Uniform temperature at varying density is not exactly representable: its false
B_t RMS converges 0.0416574, 0.000921031, 1.58422e-5, 2.53764e-7 T/s. Restoring
the full rounded E in the finest case gives a Faraday replay difference of
0.00242218 relative in the physical magnetic norm, or 6.15e-10 T/s RMS. Its raw
coefficient-space relative difference is 0.0203616. These distinct conditioning
limits remain explicit despite small original-equation backward errors.

All nine transport fixtures last one complete axis transit, 1.20801 microseconds,
with prescribed density/current and a sheared analytic entropy wave. Full-domain
entropy-signal error is 0.0846755, 0.00203068, 9.31941e-5 at nr=8/16/32; physical
first-annulus error decreases to 9.65e-9. Half-dt changes the nr16 final field
by only 1.73e-8 of the signal. Physical entropy drift is at most 2.22e-16;
proxy/internal-energy drift and initial projection error remain separate. The
finest signed proxy/internal drifts are 3.61e-12/-1.47e-12. This is a prescribed-
current transport verification, not a force-balanced or full moving plasma.

Source b0dba42d2dd23742246776328413c3914ecfebaadf158e72d9b157c10d6485e7 and
complete records are saved in outputs/log_thermal_verification and
outputs/log_entropy_evolution. Original variable-current fc106dc1, compatible
linear 8cfb1a2f and smooth nonlinear dd060015 identities remain unchanged.
See LOG_THERMAL_NUMERICS/VERIFICATION and the hashed checkpoint for full limits.

NEXT derive the common moving state and total-energy contract before adoption:
shared conservative continuity and bulk entropy advection, true conservative
velocity/vorticity, momentum/drift stress, density/pressure work and complete
current/enthalpy/Poynting transport. Resolve the two-nodal-entropy versus old
cell-energy partition explicitly; then verify independent parent-species rates
and full nonlinear wave/contact trajectories. Keep frc-7mm, frc-accuracy and the
overall goal open. Heated free edge/exact vacuum, physical transport/ionization,
real ends/coils, shocks, kinetic validity, converged formation and UI remain.
Work and issue state are local; there is no configured remote or publication.
