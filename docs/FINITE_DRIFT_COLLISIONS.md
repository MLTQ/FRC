# Paired finite-drift collision progress

## Moving paired closure integrated into the main checkout (2026-10-16)

After original run 58668 exited successfully, the exercised paired-collision
solver, experiment configs and existing field/physical-comparison tools were
integrated from `codex/finite-drift-collisions` into the main checkout. This
preserves the earlier complete trajectories and all newer atomic/radiation
work. No live source was changed. The resulting experiment source identity is
the same `31721f1fef01b747862f5263f3ddb92e5dfbca0d3fc1995e5e03ebf229d17afc`
as the archived paired runs. Prescribed-resistance mode remains explicit and
available; paired mode still replaces both prescribed resistance and separate
temperature exchange with one force/heat law, without double counting.

Paired 16000-step run 16339 remains live in `/private/tmp/frc-finite-drift`;
its latest saved diagnostic at this checkpoint was 12800/16000, not a completed
pulse. Poll that exact handle and use the now-integrated
`examples.material_coil_plot --configured-fields` for the full 8000/16000
comparison after successful termination. Keep this private worktree unchanged.

Main now also runs the original-resistance strong pulse with the existing split
entropy/current integrator and numerical bulk pressure (95749); see the new
top section of MATERIAL_COMPRESSION.md. Do not change its source while live.
This is numerical-method work; it does not resolve the paired fine-grid wall
contact, physical heat flux, ionization, kinetic validity or confinement.

## Finer mesh reaches the wall-contact model limit (2026-10-16)

Spatial run 49353 terminated with exit 1 after 9473.880 s. It accepted 6147 of
8000 steps; attempted step 6148 stopped with: "Plasma reached tube:
wall/sheath/contact physics is not implemented". The recorded source remains
31721f1fef01b747862f5263f3ddb92e5dfbca0d3fc1995e5e03ebf229d17afc.
All accepted states are retained in
`outputs/material_tube_independent_20mt_maxwellian_r6_z24_t8000`.
This is a terminal partial trajectory, not a successful full pulse and not a
reason to restart. No field, density floor or artificial continuation was added.

The last accepted time is 2.6893125 us. Its actual saved coordinates give maximum
radius 74.994630 mm in the 75 mm tube, midplane radius 42.231923 mm, maximum
half-length 388.477905 mm and particle RMS radius 19.321212 mm. Mean electron/
ion temperatures are 12.124309 / 3.101361 eV. This endpoint is reconstructed
from `step_06147.npz`, SHA256
`c866a857afd1625bc0f3afdab0a1b40bb7756c82ccb94b3c452f8a0c6db3169a`.
The summary's last periodic diagnostic is step 6100, not the final accepted
state; its accumulated heat must not be labeled as a step-6147 measurement.

Existing physical sampling at three exactly common saved steps was run without
changing either solver (analysis session 8559 exited 0). Its full component
differences and state identities are in
`outputs/material_tube_maxwellian_spatial_contact_comparison/comparison.json`.
Ten-point quadrature on the union reference mesh gives 26400 samples and uses
the finer run's actual volume. This comparison ends at the accepted prefix.

| Time | B vector difference | J vector difference | Bulk velocity difference | Maximum material-position mismatch |
|---|---:|---:|---:|---:|
| 2.1875 us | 21.9179% | 96.7556% | 57.8872% | 20.7271 mm |
| 2.625 us | 62.1559% | 97.6467% | 36.4778% | 37.6754 mm |
| 2.6893125 us | 69.8034% | 93.9034% | 33.6776% | 39.6623 mm |

At the last shared time, the 4x16 run's maximum radius is only 39.851850 mm,
although its particle RMS radius is a similar 19.683057 mm. Bulk radius and
small energy residuals conceal large edge and field differences. The finer
run's last periodic work/E0 residual is bounded in magnitude by 9.9961e-5;
that does not establish the accuracy of its wall approach. Spatial and time
errors remain mixed. The actual contact time and laboratory wall impact are
not resolved predictions. Neither grid establishes FRC formation or confinement.

Collect the already-running 4x16/16000-step paired trajectory 16339 and original
trajectory 58668. Do not merge moving sources until the original is terminal.
Use these results to address edge/current/entropy spatial treatment and physical
heat/particle boundaries; simply imposing a reflective wall would introduce
new physics without establishing an applicable sheath/contact model. Keep
`frc-transport.3/.4`, `frc-vacuum.3.3` and the full goal open.

## Completed paired timestep comparison (2026-10-16)

Run33952 exited0 after all8000steps/3.5us in5172.282s. Its source remains
31721f1fef01b747862f5263f3ddb92e5dfbca0d3fc1995e5e03ebf229d17afc;
all8000 raw states and the final matching digest are present in
outputs/material_tube_independent_20mt_maxwellian_r4_z16_t8000.
Comparison13672 exited0 using the existing actual-basis comparison; its PNG
was inspected and PDF saved in outputs/material_tube_maxwellian_timestep_comparison.
Only dt changes, from0.875 to0.4375ns. The same81 common times and6400
material quadrature samples measure each physical component using actual
second-run volume weights. No temporal interpolation or state repair is used.

| Quantity | 4000 steps | 8000 steps |
|---|---:|---:|
| Final collision heat |109.01233mJ|107.87175mJ|
| Final numerical ion heat |32.00428mJ|32.33839mJ|
| Final physical ion collision heat |0.433473mJ|0.434228mJ|
| Final particle RMS radius |22.97734mm|22.95181mm|
| Final middle / maximum radius |55.26742 /62.71483mm|54.71491 /64.31159mm|
| Final mean Te / Ti |10.53349 /1.99530eV|10.45687 /2.00251eV|
| Max observed absolute external-work residual /E0 |6.12049e-4|8.95750e-5|

The small bulk changes do not establish field accuracy. Final relative vector
differences are B10.1933%, induced B10.0786%, J20.9077%, velocity2.4197%;
maximum physical-position mismatch is4.32214mm. Across the sampled history,
B reaches20.3696% at2.84375us, J62.3226% at2.58125us, and velocity18.3392%
at1.53125us. The fine run is a comparison reference, not presumed truth.
Numerical ion heating remains much larger than physical ion collision heating.

This result triggered a third full paired pulse at0.21875ns/16000steps,
live handle16339, outputs/material_tube_independent_20mt_maxwellian_r4_z16_t16000.
Its config is saved in branch60f76e4e4 at
configs/material_tube_independent_20mt_maxwellian_t16000.json. All physical
inputs, numerical pressure coefficients and scientific source are unchanged.
The separate6x24 spatial control49353 remains live. Collect complete results
and compare8000/16000 with the existing tool; do not infer convergence from
smaller energy error alone. Keepfrc-transport.3 and the full goal open.

## First complete moving paired-collision pulse (2026-10-16)

The full 4000-step pulse completed all 3.5 us with 0.875 ns steps, exact
exterior vacuum, and the declared fixed lnLambda=10 paired Maxwellian law.
Run52937 exited0 in2621.768s; all4000 raw states and the matching final digest
are present in outputs/material_tube_independent_20mt_maxwellian_r4_z16_t4000.
Its unchanged solver source is
31721f1fef01b747862f5263f3ddb92e5dfbca0d3fc1995e5e03ebf229d17afc.

The existing physical comparison now accepts different collision laws while
requiring the same initial plasma, coils, mesh, timestep and numerical pressure.
It names each changed physical input and weights drift content by actual
J.E_collision, retaining J-squared separately. Run38724 exited0; its PNG was
inspected and PDF saved in outputs/material_tube_array_collision_model_comparison.
Code/config are committed locally as f5a7bb66d on codex/finite-drift-collisions.

| Quantity at3.5us | Prescribed eta=1e-4 ohm m | Paired Maxwellian lnLambda=10 |
|---|---:|---:|
| Particle RMS radius |25.3185mm|22.9773mm|
| Midplane surface radius |61.6000mm|55.2674mm|
| Maximum surface radius |64.7229mm|62.7148mm|
| Mean Te / Ti |13.7391 /1.0639eV|10.5335 /1.9953eV|
| Total collision heat |238.0457mJ|109.0123mJ|
| Numerical ion heat |15.6337mJ|32.0043mJ|
| Max observed absolute work residual /E0 |4.5985e-5|6.1205e-4|

The paired final RMS radius is9.247% smaller, but both exceed the initial
18.7082mm. Paired minimum recorded midplane radius is32.7585mm. This is a
changed exploratory trajectory, not convergence, formation or confinement.
Paired electron/ion collision heat is108.5789 /0.43347mJ: the separate32.0043mJ
numerical ion heat is substantial and cannot be presented as physical collision
heating. The prior63.449mJ hypothetical drag integral on the old prescribed-eta
trajectory did not predict the new109.0123mJ coupled result.

At peak drive1.26875us, independent exterior/volume linkage differences are
+0.2057% for bias and-0.3863% for the array. Drift>=electron thermal speed
contains1.261% of particles and4.129% of actual collision power. The extreme
edge drift ratio alone does not quantify the bulk contribution.

Completed diagnostic90288 exited0; both PNGs inspected and PDFs saved in
outputs/material_array_maxwellian_collision_regime_content, diagnostic source
a7d744618c104fb346c82d3934ef04c6e7cf7a19e657c94745087710da34faff.
It explicitly labels paired power as the actual implemented law. Sparse-state
drag integral109.0713mJ differs by0.0541% from the every-step109.0123mJ observation.
This is sampling difference, not an energy repair or trajectory-accuracy metric.

Transport limitations persist: at peak drive, strong-field electron parallel
density Kn>=0.1 contains47.420% of particles and81.247% of actual drag power;
ion density Kn>=0.1 contains34.068% and76.978%. At final time the electron
parallel region contains89.199% of particles and85.016% of drag power. Over
the sampled trajectory those electron/ion regions contain70.730% /65.960% of
the drag integral. Regions overlap; conditional thermal coefficients, density
gradients and fixed Maxwellians do not establish kinetic validity or heat loss.
The maximum local collision relaxation rate times dt is0.07221, only one of
the system's timescales. Actual spatial heat transport remains frc-transport.4.

Paired time control33952 (8000steps,0.4375ns) and new spatial control49353
(6x24 plasma,6 vacuum cells/12modes,8000steps) remain live on the isolated
worktree. Full time/spatial comparisons are pending. Do not change their
scientific source. Keep frc-transport.3 and the overarching goal open.

## Actual density scales expose bulk transport limits (2026-10-16)

The existing collision diagnostic now evaluates every actual plasma quadrature
sample of the complete original prescribed-eta array trajectory. It combines
conditional [NRL thermal coefficients](https://www.govinfo.gov/content/pkg/GOVPUB-D210-PURL-LPS118567/pdf/GOVPUB-D210-PURL-LPS118567.pdf)
(printed pp.34,37–39) with the physical density gradient and actual B direction.
These near-thermal, weak-coupling estimates do not supply a finite-drift cutoff
or guarantee a local transport closure. Density scales omit temperature/field
gradients, cell-interface jumps and nonlocal distributions.

For the actual Q1 meridional map A=dx/da, det(F)=(R/r)det(A), so
`grad_x ln(n)=A^-T[grad_a ln(n0)-grad_a ln(det(F))]`. The configured initial
profile gives grad ln(n0) exactly; analytic bilinear-map mixed derivatives
supply the remaining terms. No cloud interpolation or density floor is used.
The initial gradient agrees with its independent closed form to2.02e-16
relative L2. Independent central differences of saved density/coordinate maps
at initial, peak and final states reduce error from2.76–2.79e-6 to6.91–6.98e-7
when the differencing interval halves. This supports the diagnostic derivative,
not trajectory accuracy or model validity. Formula/recipe and values are saved.

Review masks use existing conservative thresholds: thermal log<5, or mean free
path divided by density scale >=0.1. Electron weak/strong-field groups separate
at omega_ce/nu_ei=1; in the strong group the density derivative is along B.
Regions overlap. Ratios are warnings to examine physics, not validity theorems.

| Actual saved state | Review region | Particle fraction | Candidate drag-power fraction |
|---|---|---:|---:|
| Initial | Weak-field electron density Kn>=0.1 | 99.179% | undefined, zero current |
| 1.26875 us | Strong-field electron parallel density Kn>=0.1 | 47.120% | 78.861% |
| 1.26875 us | Ion isotropic density Kn>=0.1 | 30.794% | 67.665% |
| 3.5 us | Strong-field electron parallel density Kn>=0.1 | 92.720% | 93.605% |

Initial unavailable thermal coefficients occupy only3.25e-9 of particles;
small positive thermal logs occupy1.553%. Thus the observed transport-scale
problem extends far beyond the dilute cold edge. Along the unchanged original
trajectory, strong-field parallel review regions account for61.994% of the
hypothetical fixed-log10 drag integral, ion density review regions40.398%,
and small thermal logs2.923%. These are path sensitivities, not predicted heat
losses or a new coupled trajectory. Candidate power differs from original eta J².

Full diagnostic10579 and independent derivative calculation83409 exited0.
Evidence: outputs/material_array_collision_regime_content; source
68546fde89b65e5a2e705aba9ba766b281888318c794c9b841420384a76998da.
The new scientific region-content PNG was inspected; PDFs generated. Code is
local branch commit50bdadad0, with running solver source31721f1f...d17afc unchanged.
No unit suite or new verification framework. Newfrc-transport.4 tracks physically
derived directional/nonlocal heat evolution and explicit vacuum-edge/end energy
assumptions. Paired drag/species exchange has no spatial heat flux and cannot
resolve this by itself. Complete paired-pulse/time comparisons remain first;
repeat the region analysis on their actual states before choosing a heat closure.

The current formation exploration combines prescribed scalar resistance with
small-drift electron-ion temperature exchange. In the completed long-array run,
7.985% of instantaneous peak-drive Joule power lies at drift above the electron
thermal speed. A physically paired drag/heat law is therefore a concrete next
simulation improvement, even if the separate numerical refinements converge.

`plasma_bench/maxwellian_pair_collisions.md` derives the implemented force and
both thermal powers from binary-scattering moments of two displaced isotropic
Gaussians. The kinetic reference is
[Pfefferlé, Hirvijoki & Lingam (2017)](https://arxiv.org/abs/1701.08037),
especially the Landau operator, Gaussian potentials and momentum/energy identities.
The implemented specialization supplies instantaneous Maxwellian moments, not the
paper's complete Hermite system or a conductivity tensor. Coulomb logarithm is an
explicit positive input. No physical cutoff, floor, fitted drag or atmosphere is
introduced. Heat can transfer out of either species; it is not clipped.

## Complete collision-source trajectories

`outputs/maxwellian_pair_relaxation` contains six full homogeneous relaxation
runs: equal masses and physical electron/argon masses, each at 100/200/400 RK4
steps over 100 initial small-drift relaxation times. Both three-component species
velocities and temperatures evolve. Initial relative speed is three combined
thermal speeds, with a nonzero common boost and declared lnLambda=10.

Independent reference integrates binary momentum and energy exchange directly
over relative velocity space, using the conditional Gaussian pair centre velocity
for the other three dimensions. It uses no production moments, erf or incomplete
gamma. Conservation alone predicts eventual common velocity and temperature.
No magnetic field, spatial transport or plasma-vacuum geometry belongs to this
collision-only source problem; it is not a full formation trajectory.

| Case | Steps | Max relative drift error | Max temperature error | Max energy drift / E0 |
|---|---:|---:|---:|---:|
| Equal mass | 100 | 3.6541e-5 | 9.8291e-5 eV | 1.9656e-5 |
| Equal mass | 200 | 2.2497e-6 | 5.8310e-6 eV | 1.0909e-6 |
| Equal mass | 400 | 1.3895e-7 | 3.7014e-7 eV | 6.4183e-8 |
| Electron/Ar | 100 | 3.0163e-5 | 8.9363e-4 eV | 1.8588e-5 |
| Electron/Ar | 200 | 1.9783e-6 | 6.4300e-5 eV | 1.0234e-6 |
| Electron/Ar | 400 | 1.2853e-7 | 4.4149e-6 eV | 5.9641e-8 |

The error scales approximately fourth order. Fine momentum errors normalized by
initial drift momentum are below 3.2e-15, and entropy is nondecreasing in these
saved histories. Equal masses end at 5.000000321 eV against the independent 5 eV
equilibrium. Physical electron/argon temperatures end at 39.98627/1.01336 eV over
1.08804 us, still far from their eventual common 20.49981 eV. Current relaxation
and temperature equilibration have different timescales. No energy reset is used.

Run 53955 exited 0; frozen source
79b7bf4c446736318c5f970ca50bcce508a9bf269ed27bdd02477f3d4e889080.
All accepted states and independent reference values are saved. PNG inspected;
PDF generated. No unit-test suite or recursive verification framework was used.

## Sensitivity on unchanged pulse states

`outputs/material_array_collision_sensitivity` evaluates the candidate on the
completed long-array run, without modifying or evolving any state. The declared
lnLambda=10 is an exploration, not a justified physical screening prescription.
Integrals along the original saved path are 236.381 mJ for prescribed eta,
77.678 mJ for the small-drift Maxwellian estimate and 63.449 mJ with finite drift.
The latter is 18.3% below its corresponding small-drift estimate. These values
are **not predictions of heating in a newly evolved pulse**. The original solver's
every-step Joule integral remains 238.046 mJ; the smaller 236.381 value here uses
the much sparser diagnostic-state quadrature.

At recorded peak drive 1.26875 us, powers are 190.270 kW prescribed and 33.329 kW
candidate. Maximum local relative-relaxation rate times the original 0.875 ns
step is 0.07222 across sampled times; this screens one collision timescale only.
It does not establish stability of the coupled equations. The candidate's actual
electron/ion powers sum to the lost drift power without a thermal correction.
Run 61137 exited 0, source
ae11182af08eb9de9989fd250a5113268ee9ca50d4c75396e7d707fdc53c6810.
The scientific PNG was inspected; PDF generated.

## Moving-plasma integration now running

`frc-transport.3` remains in progress. Commit `279bf0a8f` on local branch
`codex/finite-drift-collisions` implements one paired collision hook in the sparse
worktree `/private/tmp/frc-finite-drift`. It is not yet merged into the main
checkout: three original runs still require unchanged main scientific source.
All new physical output is saved under `/Users/max/Code/FRC/outputs`, not in the
temporary checkout. The local Git branch preserves the implementation.

`material_collision_sources` returns R_e/(e*n) to MaterialVectorCurrent, whose
driven child reuses exactly that forcing. The same call supplies Qe/Qi through
the existing dual entropy work pairing, applied once. Paired mode requires zero
prescribed eta and disables separate thermal exchange. Full physical current,
mass, Hall, stresses, fluid, pressure, exact vacuum and source work remain.
Signed species heat integrals, collision entropy and local relaxation rate are
observed separately from numerical ion heat and the full energy/work residual.

Two complete 3.5 us array trajectories are requested, both with declared fixed
lnLambda=10, the existing 4x16 mesh and numerical bulk viscosity C1=0.5, C2=2:

- LIVE handle 52937: `outputs/material_tube_independent_20mt_maxwellian_r4_z16_t4000`,
  0.875 ns integration, recorded 150/4000 steps at 102.49 s elapsed.
- LIVE handle 33952: `outputs/material_tube_independent_20mt_maxwellian_r4_z16_t8000`,
  0.4375 ns integration, first accepted step recorded.

Both actual handles were polled live and use frozen source
`31721f1fef01b747862f5263f3ddb92e5dfbca0d3fc1995e5e03ebf229d17afc`.
The early 0.13125 us coarse state has electron/ion collision heat
-21.7320458/+21.7320466 uJ, predominantly exchange before the pulse begins.
External-work residual/E0 is 6.67e-15 there. This is early integration evidence,
not a completed pulse, accuracy claim or prediction of formation. All accepted
states and failures remain saved. No unit tests were added or run.

Next collect the exact live handles; do not restart on observation timeouts or
edit their worktree source while running. Compare complete pulse trajectories
and conservation with the half-timestep control and original prescribed-eta
run. Spatial sensitivity remains required. Integrate the code into main only
after all three original-source physical jobs have terminal evidence.

Even exact instantaneous Maxwellian moments do not guarantee distributions stay
Maxwellian. Stress anisotropy, heat flux, kinetic instabilities, magnetized
collision corrections, screening, neutrals, ionization and preparation remain
unresolved. Full moving trajectories and their numerical sensitivity are pending.
