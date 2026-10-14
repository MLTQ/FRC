# Evolving the current at the driven boundary

The missing wall-current component is now part of the physical state. Two
complete 24 µs pulse runs finish. On the same eight-cell grid, the space-time
current error falls from 28.06% to 0.01512% of the independent reference RMS,
and the thermal-energy discrepancy changes from -7.87% to +0.002600%.
The finer spatial calculation is still running; this is progress on the driven
component, not completion of the laboratory FRC goal.

## Physical change

The former azimuthal current basis forced Jtheta(R)=0, although the prescribed
loop voltage requires a nonzero boundary current. Add chi=phi_wall-P0 phi_wall,
where P0 is the physical L2 projection into the original interior current space.
Then chi(R)=1 and chi is orthogonal to the interior electric functions. Evolve

    J = J0(B) + chi y(z).

The wall current y has its physical mu/(n e²) inertia. It is not replaced by an
instantaneous E/eta value. Lorentz force, current advection/stress, species
transport, Joule heat, current kinetic energy and source work all use the same
full J. No compensating heat, floor, atmosphere or numerical diffusion is added.
All old implementations and their recorded sources remain unchanged.

The full Ohm equation separates into interior and wall blocks:

    M0 E0 = I00 J0_t + C y_t + f0,
    0 = C.T J0_t + D y_t + fchi.

Eliminating y_t gives Ieff=I00-C D^-1 C.T and feff=f0-C D^-1 fchi. The original
strict induction solver acts on this reduced system. Its full current rate is
then reconstructed and the complete unreduced Ohm equation checked at the
unchanged tolerance. Actual nonuniform density enters C and D at every stage.
The time integrator evolves y together with every fluid and magnetic variable.
See the wall_current_geometry/model/plasma/gauss companion documents.

## Complete end-to-end evidence

Physical conditions and independent Bessel baseline are identical to
[DRIVEN_FLUX_RESULTS.md](DRIVEN_FLUX_RESULTS.md): R35mm,L240mm,Ar+,n1e21/m³,
initial Te/Ti5/1eV, prescribed mean-field amplitude0.1mT,12µs sin² pulse plus
12µs relaxation, fixed eta38.3182µΩ m. The magnetic field, current, motion and
heating are evolved in the full moving model. The reference is linear and
constant-eta; it omits the nonlinear fluid response retained in the simulation.

| Complete run | Old current space,160 steps | Wall current,160 steps | Wall current,320 steps |
| --- | ---: | ---: | ---: |
| Space-time J RMS error / reference RMS | 28.06% | 0.01512% | 0.02048% |
| Final thermal-energy gain | 9.250815µJ | 10.041496µJ | 10.041236µJ |
| Difference from reference Joule heat | -7.8717% | +0.002600% | +0.0000073% |
| Maximum Bz RMS error / peak drive | 2.79186e-4 | 2.79185e-4 | 2.80319e-4 |
| Maximum energy-minus-source / peak field energy | 7.56468e-5 | 7.56530e-5 | 4.73071e-6 |
| Maximum discrete/wall power discrepancy | 0.131842W | 0.000372219W | 0.000396497W |
| Maximum wall-current absolute error | imposed zero trace | 5.53866A/m² | 5.01795A/m² |

Continuum Joule heat is 10.041234979µJ. Maximum particle-number change in the
new runs is below7.8e-17 relative. The maximum instantaneous total-power minus
source discrepancy is below3.2e-12W. Energy drift is scaled by3.675µJ of peak
harmonic-field energy, not the1.33J background thermal energy.

Time refinement reduces the energy error15.9919-fold. It does not establish
current time convergence: the320-step history samples an additional early
boundary response and has a larger recorded peak current error,6.91 versus
2.25A/m², and a larger time-integrated reference error. Both use endpoint
trapezoidal time integration. The B/current spatial error remains to be measured
on the finer grid. The very small final heating difference can include
cancellation between temporal, spatial and finite-amplitude effects; it is not
a general physical error estimate or experimental validation.

Actual sessions5146 and51847 exited0 after all160/320 steps. Elapsed times were
288.00 and178.64s; peak RSS0.325 and0.319GB. Their frozen sources, accepted states,
wall-current coefficients, full stages and solver histories are saved in
`outputs/wall_current_fixed_r8_t160` and `outputs/wall_current_fixed_r8_t320`.
Direct measured comparisons are in `outputs/wall_current_comparison`.
No unit tests, automatic replay campaign or new verification framework.

## Live continuation

The full16-cell,320-step counterpart is running in actual session14381,
`outputs/wall_current_fixed_r16_t320`. Poll that handle; do not restart it or
modify its source closure. After its terminal result, compare the same physical
field/current/heating quantities, then evolve the Spitzer-resistance pulse
with the corrected boundary. Keep frc-goal.1 in progress until those results
establish the bounded component. Do not substitute another reporting milestone.

The plasma still fills a rigid periodic cylinder. There is no solved external
vacuum, physical free edge, finite-coil or real-end coupling. Weak-field Z1
transport still lacks conduction, thermal forces, viscosity, neutrals and
ionization. Neither field reversal nor axisymmetric closed contours would
prove three-dimensional stability or particle confinement. The full goal and
the separate unresolved long-contact accuracy issue remain open.
