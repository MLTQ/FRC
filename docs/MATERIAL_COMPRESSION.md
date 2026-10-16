# Numerical compression treatment and the long-array pulse

## Third original timestep: field convergence still absent (2026-10-16)

The original prescribed-resistance pulse completed 16000 steps over 3.5 us
(session 58668, exit 0, 10186.854 s). All 16000 accepted states and the matching
final digest are retained in
`outputs/material_tube_independent_20mt_viscous_r4_z16_t16000`.
Its source is unchanged:
`05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75`.

The existing actual-field comparison completed (54528, exit 0); its JSON and
inspected PNG/PDF are in
`outputs/material_tube_array_timestep_8000_16000_comparison`.
Only dt changes, from 0.4375 to 0.21875 ns. The same 81 common saved times and
6400 material quadrature samples use the second run's actual volume. There is
no interpolation in time, new fit or repair of either trajectory.

| Quantity | 8000 steps | 16000 steps |
|---|---:|---:|
| Final Joule heat |239.616087 mJ|237.561459 mJ|
| Final numerical ion heat |15.527827 mJ|15.574380 mJ|
| Final particle RMS radius |25.277822 mm|25.323649 mm|
| Final middle / maximum radius |60.796912 / 64.214075 mm|61.120178 / 64.880913 mm|
| Final mean Te / Ti |13.794141 / 1.062740 eV|13.714015 / 1.066903 eV|
| Max observed absolute work residual / E0 |1.690418e-5|9.817969e-6|

Final vector differences are B 7.5666%, induced B 7.4594%, J 18.6568% and bulk
velocity 3.4036%; maximum material-position mismatch is 8.06427 mm. Across the
sampled history, B reaches 20.7074% at 1.400 us, J 88.4446% at 1.70625 us and
velocity 30.8887% at 1.6625 us. The preceding 4000/8000 comparison gave final
B 8.2404%, J 20.5834% and position mismatch 2.46342 mm. These three trajectories
do not show settled temporal convergence. Better integrated work balance and
nearly unchanged bulk radii do not establish local field or edge accuracy.

Next method study: the existing implicit entropy/current-switch method now runs
the same strong pulse WITH the same declared numerical bulk pressure. The older
strong split failure lacked that compression treatment. New config is
`configs/material_tube_independent_20mt_viscous_split_t8000.json`, output
`outputs/material_tube_independent_20mt_viscous_split_r4_z16_t8000`, live handle
95749. It retains prescribed eta=1e-4 ohm m, 4x16 mesh, 0.4375 ns steps, all
physical terms, exact exterior vacuum and unchanged wall-contact rejection.
Its source is `31721f1fef01b747862f5263f3ddb92e5dfbca0d3fc1995e5e03ebf229d17afc`.
The discontinuous upwind pressure trace is a known numerical mechanism, not a
proven sole cause of these differences. Collect the full result before claiming
improvement; retain failure states if it stops. No fourth RK refinement was
launched, and no unit suite or new verification framework was added.

## Complete spatial comparison (2026-10-16)

Original spatial run 85554 exited 0: all 8000 steps / 3.5 us in 12445.786 s,
source 05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75.
Its final raw-state digest matches the recorded digest. The 6x24 plasma mesh,
6 vacuum cells and 12 vacuum modes replace 4x16 / 4 / 8; both use degree 3,
six-point simulation quadrature, 0.4375 ns RK4, C1=0.5/C2=2 and identical
physical inputs. Inventory normalization changes central density by 2.28 ppm.
This changes plasma and vacuum spaces together, not either in isolation.

Existing configured-fields comparison 5386 exited 0. It evaluates 81 common
physical times on 26400 common material quadrature points; actual volume of
the second run provides the weights. Read `outputs/material_tube_array_spatial_comparison`.
PNG inspected, PDF and complete comparison JSON saved. No unit tests were run.

| Quantity at 3.5 us | 4x16 | 6x24 |
|---|---:|---:|
| Joule heat | 239.6161 mJ | 247.2777 mJ |
| Numerical ion heat | 15.5278 mJ | 15.5193 mJ |
| Particle RMS radius | 25.2778 mm | 25.0652 mm |
| Midplane surface radius | 60.7969 mm | 69.3804 mm |
| Maximum surface radius | 64.2141 mm | 73.0831 mm |
| Mean Te / Ti | 13.7941 / 1.06274 eV | 13.3775 / 1.06039 eV |
| Max observed absolute external-work residual / E0 | 1.6904e-5 | 9.1335e-5 |

Final relative vector L2 differences are B 27.7235%, induced B 27.2925%,
J 60.3198% and velocity 14.9646%. Maximum material-position mismatch is
28.9723 mm (including axial and radial motion). Final component RMS differences
(r,theta,z) are B=(0.20201,0.25883,0.34610) mT and
J=(7.5283,45.0394,26.2721) kA/m2. Across the trajectory the vector differences
reach B 40.8138% at 2.8875 us, induced B 54.0387% at 2.275 us, J 86.4224% at
1.70625 us and velocity 48.8969% at 1.6625 us. This is not spatial convergence;
the already-measured timestep sensitivity remains superposed.

Independent peak-drive linkage calculation 66781 exited 0 and is recorded in
`peak_linkage.json`. At matched 1.26875 us, exterior-flux versus volume-current
linkage discrepancies improve from +1.1722% to -0.04775% for the bias pair and
from -0.62141% to -0.06267% for the long array. The refined vacuum representation
therefore improves this particular reciprocity measurement; it does not make
the full field or edge motion resolved, and no correction was applied.

The finer run's maximum radius is only 1.917 mm short of the 75 mm contact
limit, while differences between meshes are much larger. Do not interpret the
predicted clearance or cold dilute-edge motion as a lab recommendation. Similar
bulk heat and particle radius conceal large field and surface changes. Collect
the already-running 0.21875 ns controls before changing either live solver;
then prioritize the physical edge and spatial current/entropy discretization,
not another claim of accuracy based on integrated energy. `frc-vacuum.3.3`
and the full goal stay open. Paired-collision spatial/time runs are separate.

## Full isolated timestep comparison (2026-10-16)

The original long-array8000-step run88009 exited0 after all3.5us in5189.560s.
Its0.4375ns step is half the completed4000-step run's0.875ns; initial plasma,
coils, resistance, mesh, RK4 and C1=0.5/C2=2 are unchanged. Both use source
05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75.
All8000 raw states and the matching final digest are present.

Complete existing configured-fields comparison18763 exited0, with81 exactly
common physical times and6400 common material quadrature samples per time.
PNG inspected; PDF and JSON saved in outputs/material_tube_array_timestep_comparison.
Independent peak linkage calculation88092 also exited0 and saved its method
and source/input identities there as peak_linkage.json. No unit suite was run.

| Quantity |0.875ns|0.4375ns|
|---|---:|---:|
| Final Joule heat |238.0457mJ|239.6161mJ|
| Final numerical ion heat |15.6337mJ|15.5278mJ|
| Final particle RMS radius |25.3185mm|25.2778mm|
| Final midplane surface radius |61.6000mm|60.7969mm|
| Final mean Te /Ti |13.7391 /1.0639eV|13.7941 /1.0627eV|
| Max observed absolute external-work residual /E0 |4.5985e-5|1.6904e-5|

Despite similar integrated heat and bulk radius, actual final field differences
are B8.2404%, induced B8.1240%, J20.5834% and bulk velocity3.0300% in volume-
weighted relative vector L2. Final component RMS differences (r,theta,z) are
B=(5.0392e-5,1.2006e-4,9.3474e-5)T and
J=(4501.49,13995.31,14843.34)A/m2. Maximum material-position mismatch is2.4634mm.
Across sampled times, vector B differences reach14.5805% and J71.9883%, both
at1.575us. Induced B reaches35.4884% at2.275us and velocity21.9592% at1.6625us.
Early relative component differences near zero amplitude are ill-conditioned;
the saved absolute differences remain necessary. These are common material-
label comparisons; each state's physical coordinates differ.

At matched peak drive1.26875us, bias exterior/volume linkage discrepancy changes
from+1.2567% to+1.1722%, array from-0.6184% to-0.6214%. Timestep refinement does
not resolve the finite-space reciprocity discrepancy. Small energy residuals
also do not certify detailed field accuracy. Two timestep levels do not establish
convergence; the visible current structure remains materially time sensitive.

Consequent next numerical action is a full third timestep level: new config
configs/material_tube_independent_20mt_viscous_r4_z16_t16000.json changes only
the timestep to0.21875ns (title descriptive), retaining the full3.5us pulse.
Actual handle58668 runs it in the main checkout; source remains05555dc9...b3c75.
The first attempted launch exited before simulation because the assumed config
filename was absent; the input was then copied from the completed8000-step
run's actual config. No simulation was restarted or overwritten.
Original spatial refinement85554 remains live; compare it against the completed
8000-step4x16 run when terminal. The paired-collision time/space studies have
their own unchanged source in the isolated worktree. Keep frc-vacuum.3.3 open.

## Complete matching bias control (2026-10-16)

Bias-only run 46579 exited 0 after all 4000 steps / 3.5 us in 2539.966 s.
It has the same initial plasma, 0.875 ns timestep, C1=0.5/C2=2 numerical
compression pressure and source identity 05555dc9...b3c75 as the completed array.
All 4000 accepted raw states remain; the final saved-state digest matches.
The existing comparison 51962 exited 0; PNG inspected and PDF generated in
outputs/material_tube_array_vs_bias_viscous_comparison.

| Quantity at 3.5 us | Bias only | Nine-loop 20 mT pulse |
|---|---:|---:|
| Particle RMS radius | 27.4959 mm | 25.3185 mm |
| Midplane surface radius | 54.2617 mm | 61.6000 mm |
| Maximum surface radius | 54.2617 mm | 64.7229 mm |
| Mean Te / Ti | 3.4846 / 0.3756 eV | 13.7391 / 1.0639 eV |
| Joule heat | 0.0005202 mJ | 238.0457 mJ |
| Numerical ion heat | 1.95e-13 mJ | 15.6337 mJ |
| Max absolute external-work residual / E0 | 2.1283e-11 | 4.5985e-5 |

The array's final particle RMS radius is 7.919% smaller than bias-only, but both
exceed the initial 18.7082 mm. The array's midplane surface briefly reaches
33.9148 mm versus its initial 35 mm, then expands beyond the bias-only surface.
The minimum sampled particle RMS remains its initial value for both runs. A
local surface pinch is not bulk confinement. Reversed axial fields and plotted
flux contours do not establish resolved FRC topology, stability or confinement.
The bias field is constant, so the comparison's peak-source linkage at t=0 is
identically zero and supplies no nontrivial reciprocity evidence for that run.
The array's independent peak linkage differences remain +1.257% / -0.618%.

## Complete two-loop numerical sensitivity (2026-10-16)

The older successful 1.75 ns two-loop pulse without numerical compression
pressure is compared with its completed 0.875 ns C1=0.5/C2=2 counterpart.
Both numerical choices changed; this is not isolated time convergence. The
extended existing examples.material_coil_plot tool evaluates actual bases at
81 common saved times with ten-point Gauss quadrature on the union material
mesh. All components and physical-position differences are retained. Its code
is committed as 350f14014 on codex/finite-drift-collisions, isolated from live
main-checkout sources; the actual experiment solver source is unchanged.

At final 3.5 us, relative vector L2 differences (weighted by the second run's
actual volume) are total B 14.958%, induced B 14.696%, J 25.224% and bulk velocity
10.682%. Maximum material-position mismatch is 4.458 mm. Across the sampled
history, vector differences reach B 19.901%, J 40.160% and velocity 55.301%.
Relative component differences near zero fields can be large; absolute component
RMS differences are also saved. Small energy errors, 2.7201e-5 and 1.1540e-5,
do not establish local field accuracy. The running fixed-treatment time and
space controls remain necessary before formation interpretation.

Completed headless comparison 92245 exited 0; scientific PNG inspected, PDF
saved in outputs/material_tube_two_loop_numerical_sensitivity_agg. The first
postprocessing attempt 96694 exited 134 without results; its failure record is
preserved in outputs/material_tube_two_loop_numerical_sensitivity. A native
rendering problem was suspected, not proven. Simulation trajectories were
unaffected. The tool now saves measurements before plotting. No unit tests,
exhaustive replay or new verification framework were added.

This is an optional solver change for strong compressive flow, not a physical
transport coefficient or an established FRC formation model. The original
nonviscous trajectories remain unchanged and preserved. The option defaults off.

## What the actual strong pulse exposed

The independent nine-loop20mT pulse failed at1.064us with1.75ns RK4 steps;
its last accepted state was already numerically explosive. The split-current
method also failed after613steps. A complete restart with0.875ns steps instead
reaches1.101625us before a stage folds a material cell. A short0.4375ns continuation
from the earlier step607 reaches a similarly skewed cell near an end. It does
not establish full time convergence because it shares the earlier coarse prefix.

Saved original stages show field/current growth before the entropy explosion;
Coulomb temperature exchange is much slower on that step. The restricted
fixed-fluid EM eigenvalues do not demonstrate a simple pure electron-gyro
stability limit. Full-system stability and spatial accuracy remain unresolved.
Read outputs/material_independent_pulse_failure_analysis and the complete failed
outputs/material_tube_independent_20mt_r4_z16_t4000; no states were repaired.

## Explicitly numerical pressure and heat

For theta=div(v), c²=sum(gamma_s*p_s)/rho and minimum physical cell altitude h:

    u = h max(-theta,0)
    q = rho (0.5 c u + 2 u²)
    force = integral q div(delta x) dV
    heating = -q theta >= 0

The same cylindrical geometry/quadrature pairs the force and heat. The force
changes centre-of-mass momentum; the declared heat partition puts the loss in
ions through the actual dual entropy representation. Numerical heat is reported
separately from eta J². It does not alter physical resistance or add relative
current drag. Exterior matter/conductivity stays exactly zero. In smooth flow
q decreases as O(h)+O(h²), and q=0 in expansion. This scalar treatment cannot
resolve all shear/hourglass modes or guarantee an untangled mesh. Energy after
finite RK4 steps remains an observation, never a fitted correction.

The method uses the artificial-pressure idea of von Neumann and Richtmyer,
[1950 paper](https://doi.org/10.1063/1.1699639), and conjugate force/work reasoning
discussed in [Caramana et al.1998](https://cnls.lanl.gov/~shashkov/papers/comp.pdf).
It does not implement that paper's complete staggered algorithm or claim its
exact fully discrete conservation; the formula above is the implemented choice.

## Complete independent wall-shock baseline

The same material fluid, mass and dual heating evolve a uniform gamma5/3 gas
with Ar ion plus electron mass, n=1e18/m³ and10km/s incoming flow. Initial total
pressure is0.01 rho u0². A fixed end wall stops the downstream fluid; rigid radial
constraints give a one-dimensional shock. The moving upstream boundary's reaction
work is integrated using the actual RK stage weights. No electromagnetic or
species exchange terms belong to this Euler fixture.

Mass/momentum/energy jumps give shock speed3457.183m/s toward the inflow,
compression3.89252846 and p2/p0=135.5718293. At20us the exact front isz=0.43085634m.
These values are independent of solver output. All three complete trajectories
finish with the following volume L1 errors:

| Axial cells | dt, ns | Density | Velocity | Pressure | Energy/work residual / E0 |
|---:|---:|---:|---:|---:|---:|
| 50 | 50 | 3.1206% | 1.6184% | 10.4266% | 0.0045% |
| 100 | 25 | 1.9824% | 0.9162% | 6.1963% | 0.1049% |
| 200 | 12.5 | 0.9623% | 0.4845% | 2.9564% | 0.0523% |

The inspected PNG (also saved as PDF) shows narrowing shock width and correct downstream
pressure, with a localized density/wall-heating defect still visible. The first
50-cell mesh uses its radial altitude as h; on finer meshes h is axial. The
energy residual is not monotone from50 to100, and the given timestep refinement
alone does not remove heating-integration error. A complete100-cell run at four
times as many steps (3200steps,6.25ns) reduces the relative energy/work residual
from1.04936e-3 to3.08385e-6 (340-fold), with density/velocity/pressure L1 errors
1.8912%/0.9074%/5.7313%. It retains the correct shock profile. The inspected
PNG and generated PDF are in outputs/material_wall_shock_viscosity_t4, source
be052c85250c68fe201c8832b6e4d2e6ea3a5940b101b28a3a51dbb93745f030.
This isolates a substantial timestep contribution; the remaining wall-density
defect and shock width are not removed by that time refinement.

Evidence: outputs/material_wall_shock_viscosity_b; frozen source
112a2fde261fde404b4de8c1cd0755a37f330af0cb60acbb88bc1938df26365b.
The earlier single-radial-cell setup failed the mesh contract before evolution
and remains in outputs/material_wall_shock_viscosity. No unit tests were used.

## Physical pulse in progress

LIVE13616: outputs/material_tube_independent_20mt_viscous_r4_z16_t4000.
Same complete initial plasma, independent bias/long-array waveforms and0.875ns
step, now with C1=0.5,C2=2 numerical bulk pressure. Frozen source
05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75.
By recorded1.40us (1600/4000steps), this run has passed the previous fold time
without rejecting the mesh; it is still unfinished. Middle/max radii34.476/34.943mm,
Joule heat132.433mJ, numerical ion heat10.900mJ and work residual/E0=5.2092e-6.
Collect its actual terminal state; preserve failures and do not edit live source.
Matching two-loop control84694 is also LIVE in
outputs/material_tube_two_loop_20mt_viscous_r4_z16_t4000, same source, timestep and
numerical viscosity. Compare these to isolate geometry from numerical changes.
A completed damped trajectory would still need physical comparison and spatial
assessment. It would not establish ionization, transport validity, confinement,
three-dimensional stability or a laboratory design.

## Prescribed field extent

The independent circular-loop axis formula evaluated on these configured inputs
gives-19.9mT at centre in both arrangements. Atz=0.35m the two-loop source is
-0.63812mT and the long array is-17.14213mT:3.21% versus86.14% of centre magnitude.
Their shared100uT central bias falls to3.2066uT there. The inspected figure in
outputs/material_source_geometry_comparison shows both vacuum fields and their
Bsource²/(2mu0) scale alongside the assumed initial on-axis pressure. Induced
screening and evolving pressure/forces are not represented by those input curves.
The result explains why equal centre fields do not give equal forcing along the
tube; it does not establish a lab coil design or force-balanced plasma.

## Long-array terminal result (2026-10-16)

Run13616 completed all4000steps to3.5us in2456.917s with source05555dc9...b3c75.
Final middle/max radii61.600/64.723mm, particle RMS25.318mm, Joule heat238.046mJ,
numerical ion heat15.634mJ, meanTe/Ti13.739/1.064eV. Maximum external-work residual
relative to initial energy is4.59846e-5. The maximum final drift/thermal ratio67.87
remains large; successful time evolution does not establish transport/edge
validity, numerical convergence or formation. Matched control84694 is still live.
The new351-frame10ns-target UI uses the older completed two-loop trajectory at
the user's request; it does not yet depict or compare this new numerical case.

## Matched complete coil comparison (2026-10-16)

Control 84694 completed all 4000 steps / 3.5 us (2362.285 s), with the same source,
0.875 ns timestep, plasma inputs and numerical pressure as the long-array case.
Only prescribed coil arrangement changes. The comparison uses saved physical
states; the inspected PNG and generated PDF are in
outputs/material_tube_array_vs_two_loop_viscous_comparison.

| Quantity at 3.5 us | Two pulse coils | Nine-loop array |
|---|---:|---:|
| Midplane surface radius | 53.140 mm | 61.600 mm |
| Maximum surface radius | 67.390 mm | 64.723 mm |
| Particle RMS radius | 25.677 mm | 25.318 mm |
| Joule heat | 85.742 mJ | 238.046 mJ |
| Numerical ion heat | 2.691 mJ | 15.634 mJ |
| Mean Te / Ti | 6.996 / 0.570 eV | 13.739 / 1.064 eV |
| Maximum observed external-work residual / E0 | 1.154e-5 | 4.598e-5 |

The array produces 2.776 times Joule heat but only 1.40% smaller final particle
RMS radius. Both exceed the initial 18.708 mm RMS radius. Minimum recorded middle
radii are 34.034 mm at 1.925 us and 33.915 mm at 1.8375 us, from 35 mm initially.
The array initially suppresses the large axial variation of the outer radius,
but rebounds after the pulse. Neither surface compression nor lower RMS relative
to control establishes particle confinement. The artificial-pressure contribution
is significant: 15.634 mJ is 6.57% of the array Joule heat and 17.47% of initial
plasma energy. Its change under refinement remains to be measured.

At recorded peak drive 1.26875 us, independent exterior-flux versus volume-current
linkage differences are +0.765% for two coils, +1.257% for array bias and -0.618%
for the formation array. These finite-space discrepancies are unchanged; they do
not repair the work integral. The much smaller global work residual does not
therefore establish spatial field accuracy.

Actual quadrature content helps interpret the large drift maximum. At array peak
drive, |w|/sqrt(eTe/me)>=1 includes 0.793% of particles, 14.278% of plasma volume
and 7.985% of instantaneous Joule power. Threshold 0.1 includes 34.904% of particles
and 89.925% of power. At final time the >=1 fractions shrink to 0.0198% of particles
and 0.00173% of power. Two-loop peak >=1 fractions are 0.115% of particles and
0.847% of power. These are descriptive thresholds, not a transport-validity test;
other samples can violate collisional, isotropic or local closure assumptions.
Uniform prescribed eta makes integral J²dV weighting proportional to Joule power.
No neutral, preparation, conduction or kinetic closure has been added.

The completed array now has a separate 351-frame 10 ns-target web display at
127.0.0.1:8790/array20mt/. Its actual gaps are 9.625/10.5 ns, with no temporal
interpolation. Numerical compression coefficients and scope are disclosed.
The original two-loop /20mt/ animation is retained.

Full-pulse time and spatial refinements are LIVE 88009 and 85554, both 8000 steps
at 0.4375 ns; a matched bias-only 4000-step control is LIVE 46579. Read the newest
GOAL_STATUS boundary for paths and collection order. These runs are unfinished;
no partial prefix substitutes for completed convergence evidence.
