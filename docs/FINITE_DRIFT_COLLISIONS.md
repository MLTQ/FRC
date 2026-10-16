# Paired finite-drift collision progress

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

## Next physical integration

`frc-transport.3` remains in progress. After current physical runs finish, add one
paired collision hook to MaterialVectorCurrent and reuse its actual electrical
forcing in MaterialCoilCurrent's driven solve. Replace both old prescribed drag
and old thermal exchange for this option; do not add them twice. Species heat
must use the same dual entropy work pairing. Preserve exact vacuum and the full
current, inertia, fluid, pressure and source-work terms. Then run complete moving
trajectories with conservation and numerical sensitivity measured independently.

Even exact instantaneous Maxwellian moments do not guarantee distributions stay
Maxwellian. Stress anisotropy, heat flux, kinetic instabilities, magnetized
collision corrections, screening, neutrals, ionization and preparation remain
unresolved. The new model is not yet installed in the moving-coil simulation.
