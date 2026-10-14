# Plasma response to a prescribed axial-flux pulse

The moving collisional model now accepts a time-dependent enclosed axial flux,
including Faraday forcing and external work at both implicit time stages.
The complete fixed-resistance pulse runs finish. Magnetic-field agreement is
good, but the wall current and heating are not yet accurate enough to accept
this boundary treatment for experiment design.

## Physical calculation

The filled cylinder has radius 0.035 m, periodic length 0.24 m, fully ionized
singly charged argon at n = 1e21/m³, initial Te/Ti = 5/1 eV, and zero initial
flow, field and current. Both adiabatic indices are 5/3. Prescribe

    Bbar(t) = A sin²(pi t/T),  0 < t < T;
    Bbar(t) = 0 otherwise;    A = 0.1 mT, T = 12 µs.
    Phi(t) = pi R² Bbar(t),   E_ext,theta = -r Bbar_dot/2.

The calculation includes a further 12 µs of relaxation. The prescribed
quantity is mean axial flux, not the field at the wall. There are no finite
external coils or circuit model. The full moving finite-mass equations retain
current inertia, Hall terms, current advection/stress, species pressure,
Joule heat and Coulomb thermal exchange. No artificial diffusion, atmosphere,
floor or energy correction is added. This is a bounded mathematical case,
not a chosen laboratory geometry or an FRC formation experiment.

`AxialFluxDrive` supplies the physical inductive forcing. `DrivenModalGauss`
uses each stage's actual time and field, retaining the original residual,
boundary and positivity checks. Its source work is the Gauss integral of

    P_discrete = integral J_model E_ext dV + V Bbar Bbar_dot / mu0.

The independently reconstructed continuum wall expression is
P_wall = V Bbar_dot Bz(R)/mu0. They are distinct at finite resolution because
the old current representation forces Jtheta(R) = 0. The continuum driven
current need not vanish there. Closing the discrete energy budget does not
establish accurate wall power or Joule heating.

## Complete fixed-resistance runs

Both runs use eta = 3.8318200089006576e-5 ohm m, degree five, 16 axial cells,
half as many fluid radial elements as field cells, and 160 time steps.
All 160 states are accepted in each run. Actual sessions 4859 and 4781 exited 0.
Elapsed times are 201.76 and 836.01 s; peak RSS is 0.397 and 1.161 GB.

The independent continuum baseline solves cylindrical diffusion with physical
current inertia and prescribed enclosed flux. It uses J1 Bessel zeros and an
analytic polynomial lift that retains the nonzero wall current. See
[the equations and sources](../examples/flux_diffusion_reference.md).
The comparison is mathematical, not experimental. The reference omits the
second-order flow and temperature response, which remain in the simulations.

| Quantity over the complete 24 µs | 8 radial cells | 16 radial cells |
| --- | ---: | ---: |
| Maximum Bz RMS error / peak imposed Bbar | 2.79186e-4 | 4.82744e-5 |
| Space-time Jtheta RMS error / reference RMS | 28.06% | 21.14% |
| Final thermal-energy gain | 9.25081 µJ | 9.59284 µJ |
| Difference from continuum Joule heat | -7.87% | -4.47% |
| Maximum discrete/wall power difference | 0.131842 W | 0.0747768 W |
| Maximum energy-minus-source / peak field energy | 7.56468e-5 | 7.56475e-5 |
| Maximum speed | 0.05659 m/s | 0.05763 m/s |
| Maximum density RMS change / initial density | 8.37868e-6 | 8.65784e-6 |

Continuum Joule heat is 10.041235 µJ, integrated with 1024 modes, 256-point
radial Gauss quadrature and adaptive time quadrature across the pulse endpoint.
At selected phases, 512/1024-mode differences are at most 3.053e-13 T and
0.01812 A/m². The space-time current norm uses endpoint trapezoidal time
integration; it is not a maximum local relative error. Conservation error is
normalized by 3.675 µJ of peak harmonic-field energy, not the 1.33 J thermal
background. Maximum particle-number change is below 4.9e-17 relative.

The field improves 5.78-fold with refinement, while current improves only
1.33-fold. The heating deficit and source/wall-power discrepancy decrease,
but remain substantial. Do not infer accurate current from the magnetic
field or native energy budget. No time-refinement claim is made for this pulse.

## Evolving-resistance run

The same 8-cell, 160-step pulse with actual temperature-dependent Spitzer
resistance completes the full 24 µs. Actual session 70772 exited 0, elapsed
342.14 s, peak RSS 0.594 GB. Electron/ion mean temperatures finish at
3.02134/2.97870 eV. Resistance rises from 38.3182 to about 74.0081 µΩ m,
while the electron energy transfers mainly to ions.

The final thermal-energy gain is 7.71423 µJ, 16.61% less than the same-grid
fixed-resistance control. Maximum current at 24 µs is 0.73558 A/m² versus
52.11798 A/m² in that control. Peak speed is 0.04739 m/s. The maximum Bz RMS
control difference is 13.33% of the imposed pulse amplitude. These are
physical-model comparisons with the same unresolved wall approximation,
not numerical error estimates or a quantitative validated heating prediction.

Peak electron Hall parameter is 0.01052. Peak drift / electron and ion thermal
speeds is 7.05e-5 and 0.03093. These measured ratios satisfy the existing weak-
field/drift review thresholds; they do not establish complete transport
validity. Maximum energy-minus-source / peak field energy is 7.56462e-5;
maximum particle-number change is 5.94e-17 relative.

The [inspected response figure](../outputs/driven_flux_comparison/driven_flux_response.png)
shows the actual field delay, wall-current discrepancy, evolving resistance,
heating and motion. All three scientific workers and the plotting observer
(session 25683) are terminal, with actual exit 0 collected.

## Evidence and next physical change

Saved runs: `outputs/driven_flux_fixed_r8_t160` and
`outputs/driven_flux_fixed_r16_t160`, plus
`outputs/driven_flux_spitzer_r8_t160`. Accepted arrays, stages, solver histories
and frozen sources are retained. Direct physical comparisons are in
`outputs/driven_flux_comparison`. No unit tests or new verification framework.

Next, correct the missing driven-wall current degree of freedom together with
its finite current inertia and paired work/heating. Do not add compensating
heat or replace current only in diagnostics. Keep the original runs as
controls. The driven-plasma bead remains open.

The present rigid wall and periodic axial boundary do not provide an exact
moving plasma-vacuum interface or real ends. Weak-field Z1 transport lacks
conduction, thermal forces, viscosity, neutrals and ionization. Those issues,
finite-coil coupling and converged formation remain part of the full goal.
