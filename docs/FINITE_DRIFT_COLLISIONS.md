# Paired finite-drift collision progress

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
