# Driven-wall current: completed bounded pulse component

Four complete24µs simulations now evolve the missing wall current, including
its physical inertia, force, heating and source work. Spatial refinement reduces
full-pulse current error to0.003419% of the independent continuum reference RMS.
Time refinement reduces the energy error15.99-fold. The evolving-resistance
case also completes, with17.38% less heat than its fixed-resistance control.
This completes frc-goal.1 within its filled-cylinder scope; the full experimental
FRC goal remains open.

The [refinement and collisional-response figure](../outputs/wall_current_comparison/wall_current_refinement.png)
is visually inspected. The preceding [wall-current correction figure](../outputs/wall_current_comparison/wall_current_response.png)
retains the original/coarse comparison; its pending-spatial caption is historical.

## Physical change

The old current basis forced Jtheta(R)=0 despite a nonzero applied loop voltage.
Add chi=phi_wall-P0 phi_wall, with P0 the physical L2 projection into the old
interior theta space. Then chi(R)=1 and chi is orthogonal to interior electric
functions. The actual current is J=J0(B)+chi*y(z). The new y evolves with physical
mu/(n e²) inertia; it is not replaced by an instantaneous E/eta value.

All Lorentz, current-advection/stress, species, Joule, kinetic-energy and source-
work terms use the same full current. The full weak Ohm equations are

    M0 E0 = I00 J0_t + C y_t + f0,
    0 = C.T J0_t + D y_t + fchi.

Eliminating y_t gives Ieff=I00-C D^-1 C.T and feff=f0-C D^-1 fchi. The unchanged
strict induction solver acts on this system; reconstruction restores the full
current rate and checks the unreduced Ohm equation. Actual density enters C/D
at every stage. Gauss stages evolve all variables together at actual source
times. No numerical diffusion, floors, atmosphere, compensating heat or relaxed
acceptance bounds. Old code and recorded sources remain unchanged.

## Complete fixed-resistance simulations

All cases use R35mm,L240mm,fully ionized Ar+,n1e21/m³,initial Te/Ti5/1eV,
both gamma5/3,degree5 and16 axial cells. Prescribed mean axial field is0.1mT
sin²(pi*t/12µs), then zero during12µs relaxation. It prescribes enclosed flux,
not wall Bz or finite-coil current. Bulk motion, Hall, finite current inertia,
pressure, current transport, Coulomb exchange and electron Joule heat remain.

Eta is38.318200089µΩ m. The independent finite-inertia Bessel diffusion baseline
retains the nonzero wall current; its equations and primary sources are in
[flux_diffusion_reference.md](../examples/flux_diffusion_reference.md). The
baseline omits the nonlinear fluid/temperature response retained in the model.
This comparison is mathematical, not experimental.

| Quantity over24µs | 8 cells /160 steps | 8 cells /320 steps | 16 cells /320 steps |
| --- | ---: | ---: | ---: |
| Space-time J RMS error / reference RMS | 0.01512% | 0.02048% | 0.003419% |
| Maximum Bz RMS error / peak imposed field | 2.79185e-4 | 2.80319e-4 | 4.10337e-5 |
| Final thermal-energy gain | 10.041496µJ | 10.041236µJ | 10.041236µJ |
| Difference from continuum Joule heat | +0.002600% | +0.0000073% | +0.0000074% |
| Maximum energy-minus-source / peak field energy | 7.56530e-5 | 4.73071e-6 | 4.72980e-6 |
| Maximum discrete-source / wall-power difference | 0.000372219W | 0.000396497W | 0.0000258659W |
| Maximum instantaneous J RMS error | 2.25131A/m² | 6.91001A/m² | 1.09365A/m² |

The old8-cell/160-step model had28.06% space-time current error and7.87% too
little heat. Continuum Joule heat is10.041234979µJ. The apparent near-exact
final heat agreement can include temporal/spatial/finite-amplitude cancellation;
it is not a general physical error estimate. Energy error is normalized by
3.675µJ of peak harmonic-field energy, not the1.33J thermal background.

At the same320-step resolution, spatial refinement reduces current error5.992-
fold, peak Bz error6.831-fold and wall-power discrepancy15.329-fold. Halving dt
reduces energy error15.9919-fold. These do not imply uniform pointwise accuracy:
at the first75ns endpoint the fine current RMS error is1.26% of the then-small
reference current; it is0.250% at150ns,0.0193% at600ns and0.000163% at3µs.
Maximum wall-current error remains5.018A/m², including the unresolved fast
carrier transient. Late weak-current errors also need not improve pointwise.

The320-step history includes additional early samples. On identical150ns
samples,8-cell current error changes0.01512%→0.01428%, while its peak still
increases2.25→2.60A/m². No temporal order for current is established. The current
space/time norms use endpoint trapezoidal integration and an independent1024-
mode reference denominator. The measured spatial refinement supplies the main
current-accuracy evidence.

## Evolving-resistance physical comparison

The8-cell/320-step Spitzer pulse completes with the same corrected wall model.
Electron/ion means finish at3.021337/2.978700eV. Resistance increases from38.3182
to about74.0081µΩ m as electrons transfer energy to ions. Thermal gain is
8.295722µJ,17.38% less than the same-grid fixed-resistance case. Final maximum
current is0.73555 versus52.11797A/m²; peak speed is0.04929m/s. The largest Bz
RMS difference from the fixed-resistance reference is13.33% of the drive.
These are physical-model/control differences, not numerical errors.

Peak electron Hall parameter is0.010522; peak current-drift / electron and ion
thermal-speed ratios are7.12e-5 and0.03211. These pass existing review thresholds,
not complete transport-validity proof. Maximum energy-budget error is4.72998e-6
of peak field energy and maximum instantaneous power discrepancy2.63e-12W.
All four runs retain particle number to better than7.8e-17 relative in the
independently reconstructed delta-density integral.

## Evidence and next physical step

All scientific sessions are terminal, with actual exits0 collected:
5146 (8/160,288.00s,0.325GB peak RSS),51847 (8/320,178.64s,0.319GB),14381
(16/320,1140.54s,0.821GB),65413 (Spitzer8/320,590.66s,0.461GB). Plotting observer
97286 also exited0. Runs live in `outputs/wall_current_fixed_r8_t160`,
`outputs/wall_current_fixed_r8_t320`, `outputs/wall_current_fixed_r16_t320` and
`outputs/wall_current_spitzer_r8_t320`; comparisons are in
`outputs/wall_current_comparison`. All accepted arrays/stages and frozen
sources remain. No unit tests, automatic replay campaign or new framework.

NEXT frc-vacuum.1: outer flux must evolve consistently with a magnetic-only
exterior/source field. Preserve zero exterior fluid density/conductivity and
include exterior field energy and source work. The existing finite-volume
vacuum trace and regular-field representations require a compatible coupling.
Do not retain the closed-flux assumption that uniform-harmonic/increment cross
energy is zero when wall flux is free: for uniform B0 it is
(2*pi*B0/mu0)*integral psi_increment(R,z) dz. The fixed-footprint first step does
not complete the parent moving-interface, finite-coil and real-end requirements.

These runs still fill a rigid periodic cylinder. Weak-field Z1 transport lacks
thermal conduction/forces, viscosity, neutrals and ionization. They demonstrate
no FRC formation, three-dimensional stability or particle confinement. The
full goal and separate long-contact accuracy issue remain open.
