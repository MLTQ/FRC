# Evolving resistance, magnetic decay and plasma heating

The moving plasma model now evaluates electrical resistance from the actual
electron temperature and density at every nonlinear stage. Electron Joule heat
and electron-ion temperature exchange evolve together with the original
fluid, magnetic, Hall and finite-current-inertia equations. The implementation
is opt-in: plasma_bench/spitzer_wall_plasma.py. No artificial diffusion,
temperature floor, energy reset or changed solver tolerance was introduced.

## Completed physical runs and independent reference

The example uses a filled periodic cylinder, R=0.035m and L=0.24m, Ar+ at
n=1e21/m³, Te=5eV and Ti=1eV, both gamma=5/3. Initially v=0 and
Bz=0.1mT J0(qr), with J1(qR)=0, zero bias and fixed wall poloidal flux.
All nonlinear motion and heating remain active for12us. See the example
companion for equations, numerics and the independent scalar time integral.

The low-amplitude reference first integrates thermal equilibration, then
integrates eta(Te(t)) in the finite-current-inertia Bessel decay exponent.
It calls no production collision/RHS routine. It predicts Te/Ti of
3.288140984/2.711859016eV and eta increasing from3.831820009e-5 to
6.630526087e-5 ohm m: a73.04% increase. Final magnetic amplitude is
0.002273558842 of initial, compared with0.012474446038 if eta were held
at its initial value. The changing-resistance field isabout5.49 times smaller.

| Integrator / radial cells / steps | Final Bz L2 error / remaining signal | Final Jtheta L2 error / remaining signal | Physical energy drift / initial field+current energy |
|---|---:|---:|---:|
| Original / 8 / 80 | 3.27565e-5 | 7.73256e-6 | 4.65766e-4 |
| Original / 16 / 80 | 8.07015e-6 | 4.75221e-5 | 4.65773e-4 |
| Coordinate addition / 8 / 80 | 3.27567e-5 | 7.77363e-6 | 4.65774e-4 |
| Coordinate addition / 8 / 160 | 3.33292e-5 | 4.57165e-5 | 2.91367e-5 |

Radial refinement improves the final B error by4.06 times. Current error
does not improve; no formal current-convergence claim is made. Initial
reconstruction errors are retained in the raw evidence, including B errors
8.65601e-5/3.22284e-6 and current errors5.41466e-4/3.44981e-5 on these grids.
All comparisons use independent order20 physical reconstruction. The magnetic
reference is linear, while the simulation includes nonlinear current forces
and spatially varying Joule heat. The matching coordinate-addition80/160
runs reduce energy drift by15.9858 times. B/J errors do not time-converge in
these measurements; spatial, finite-amplitude and small numerical departures
remain. The final original-versus-coordinate-addition80-step vector B
difference is1.90935e-10 of initial Bz in independent physical L2 norm.

The initial field-plus-current energy is5.96342e-7J against1.33183J of thermal
energy. Energy drift is deliberately normalized by the small field budget;
the table is not a thermal-background-relative conservation claim. Mean
Te+Ti rises by about2.688e-6eV as magnetic energy becomes heat. Native energy
and independently reconstructed physical energy are recorded separately.

## Assumptions and evidence

This is a conditional Z1 low-field scalar Spitzer closure with coefficient
0.51, Maxwellian equilibration, and leading-me/mi electron Joule heating.
The [NRL Plasma Formulary 2023, pp.36–38](https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf)
supplies the transport convention and its restrictions. Initial maximum
electron Hall parameter is about0.00831; electron/ion mean free paths divided
by1/q are0.04849/0.00415. Collision times divided by the initial field-decay
time are0.000173/0.00891. Coulomb logarithms are modest, and these ratios do
not certify a full transport model. Hall/current inertia remain in the solver.
The scalar closure is not the strongly magnetized transverse resistance.

Conduction, thermal forces, frictional heat flux, viscosity, neutral gas and
ionization remain absent. Their omission is compatible with the first-order
uniform-temperature reference, but not a complete nonlinear laboratory model.
There is no free plasma-vacuum interface, external coil drive or physical end
model here. This does not close frc-7mm.7 or demonstrate FRC formation.

Raw states, stages, source closures and solver histories:
outputs/spitzer_bessel_decay_r8_t80 and _r16_t80. Actual sessions61936 and53684
exited0 and were collected; elapsed139.49/543.83s, peak391561216/1111425024B.
The inspected figure outputs/spitzer_bessel_comparison/temperature_dependent_decay.png
uses the completed8/80 run; its source record is beside it. No unit tests,
new verification framework or exhaustive replay was added.

Original time refinement8/160 failed at step17 after1.2us on the unchanged
wall constraint. Actual58978 exited0 with success=false and was collected.
An observational replay92257 reproduced the rejection and saved arrays/stack:
normalized trace5.87858e-14 exceeds the5.68434e-14 bound after cancellation of
two tiny constrained vectors. Absolute trace norm1.75477e-28, coefficient
norm4.14217e-16. This is preserved failure evidence, not a completed transit.

The opt-in LiftedModalGauss constructs magnetic fields after adding their
independent coordinates, with original equations and bounds. Both complete
12us8/80 and8/160 runs now pass. Actual sessions99314/35233 exited0 and were
collected; elapsed122.10/110.38s, peak384139264/367591424B. All scientific
workers are terminal. The original frozen implementation/results remain
unchanged. Use examples/lifted_spitzer_decay.py for this corrected arithmetic.

Evidence: outputs/lifted_spitzer_decay_r8_t80 and _r8_t160, plus comparison
summary.json in outputs/spitzer_bessel_comparison. The maximum particle-number
change over the four completed decay runs is5.13e-17 relative. No source or
acceptance bound was changed during any run. frc-transport.2 is complete for
this bounded component; full transport and FRC goals remain open.

Next: frc-goal.1 couples a prescribed axial-field pulse to the collisional
plasma with explicit inductive forcing and external energy input, starting
with an independently solvable cylindrical diffusion limit. It is an actual
driven simulation milestone, not another reporting/replay framework.
