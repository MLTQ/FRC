# A driven plasma coupled to a magnetic-only exterior

Five full 24 µs simulations complete with an evolving plasma boundary flux.
At the source peak of 100 µT, the plasma mean field is 51.914 µT with fixed
resistance and 57.625 µT with evolving resistance. The latter deposits 8.78%
more heat in this prescribed-source experiment. This differs from the earlier
prescribed-enclosed-flux control and demonstrates why the physical boundary
matters. These are mathematical comparison fixtures, not formation experiments.

[Completed response figure](../outputs/vacuum_pulse_comparison/vacuum_coupled_response.png)
was visually inspected. All raw accepted states, stages and frozen sources are
saved alongside the measured histories. No unit tests or new verification
framework were added or run.

## Physical implementation

The new VacuumCylinderGeometry includes outer psi and theta electric/current
values as physical degrees of freedom. Faraday evolves the boundary flux;
weak Ampere uses the derivative of interior plus exterior magnetic energy.
Current inertia, forces, advection/stress, species transport, Joule heating and
current kinetic energy all use this same current. No auxiliary wall-current
state, artificial atmosphere, exterior density/conductivity, numerical diffusion,
energy repair or relaxed nonlinear/induction tolerance is introduced.

The exterior extends to infinity in r. Its analytic K-Bessel radial elimination
uses the same polynomial axial trace as the plasma. All axial modes participate;
this is a Galerkin vacuum, whose continuum error is measured below. It is not
pointwise curl-free at finite axial resolution. The constant induced-flux mode
has zero exterior field. The source is an ideal cylindrical current sheet at
C=75 mm; its imposed field is uniform inside C and zero outside. Induced fields
extend through it, with no extra conducting wall.

Total magnetic energy includes the full induced energy, source background energy
pi C² L B_source²/(2 mu0), and the nonzero cross term
(2 pi L/mu0) B_source mean(psi_ind(R)). Source work is computed independently from
source current and linked-flux rate, including the plasma contribution. In the
uniform run the cross term is -3.5343 µJ at the source peak and cannot be omitted.
The source radius determines background energy/work, not an imposed plasma flux.

All original motion, Hall, finite carrier masses, species pressure, Coulomb
exchange and electron Joule terms remain. The existing driven Gauss algorithm
is reused with an explicit free-flux field type and coordinates. Its stage and
endpoint bounds are unchanged. The new modules have equation/contract companions.

## Completed physical comparisons

Conditions: plasma R=35 mm, periodic L=240 mm, fully ionized Ar+, n=1e21/m³,
initial Te/Ti=5/1 eV, both gamma=5/3. B_source is a 100 µT sin² pulse over 12 µs,
followed by 12 µs relaxation. Fixed eta is 38.3182 µΩ m. A separate uniform case
starts with zero field. Other cases initialize a 10 µT axial magnetic mode that
matches the continuum vacuum; this supplies nonzero exterior field and energy.
Any topology already present in that seed is initialized, not formed.

The independent fixed-resistance pulse uses J0-zero cylindrical diffusion modes
with physical current inertia. The axial seed uses the independent continuum
Robin eigenvalue and K-Bessel exterior, not production's discrete axial modes.
Equations and primary sources are in
[vacuum_pulse_reference.md](../examples/vacuum_pulse_reference.md).
The nonlinear simulation retains fluid response omitted by that linear baseline.

| Seeded fixed-resistance run | 8 radial / 16 axial, 160 steps | 16 radial / 32 axial, 160 steps | 8 radial / 16 axial, 320 steps |
| --- | ---: | ---: | ---: |
| Space-time current RMS error / reference RMS | 0.005962% | 0.001085% | 0.005512% |
| Space-time vacuum-vector RMS error / reference RMS | 0.009513% | 0.0003213% | 0.009513% |
| Maximum Bz RMS error / peak drive | 4.65181e-5 | 4.66985e-6 | 4.65220e-5 |
| Maximum reconstructed interface Bz jump / drive | 1.94465e-4 | 2.40533e-5 | 1.97036e-4 |
| Maximum energy-minus-source / peak plasma field energy | 7.56304e-5 | 7.56304e-5 | 4.73562e-6 |
| Thermal-energy gain | 2.4231033 µJ | 2.4231034 µJ | 2.4228428 µJ |

Joint radial/axial refinement reduces current error 5.50-fold, vacuum error
29.60-fold, peak Bz error 9.96-fold and the interface jump 8.08-fold. Halving
the time step reduces energy error 15.97-fold. Current/field error is mostly
spatial here; no fourth-order current convergence claim follows from the energy.
The energy scale is 3.675 µJ, not the approximately 1.33 J thermal background.

Current norms use physical volume quadrature. Vacuum-vector errors use equal-
weight samples at R, 1.5R and 2R and cellwise axial Gauss nodes; the continuum
reference harmonic has exactly the same mean-square value on these periodic
samples as its continuous axial mean. Space-time norms use endpoint trapezoids.

Global norms do not establish uniform startup accuracy: the uniform coarse run
has 4.31% current RMS error relative to its small reference current at the first
150 ns endpoint. The seeded fine run's first relative error is 0.0635%, but its
nonzero initial seed changes that denominator. Fast carrier transients are not
fully time-resolved. Reference mode doubling changes sampled J by at most
2.02e-5 A/m² and Bz by 3.69e-16 T; this does not bound all model errors.

Independent seeded continuum Joule heat is 2.4228261 µJ. The 320-step thermal
gain differs by 0.000688%, with finite-amplitude/temporal/spatial contributions
possible. Maximum instantaneous power mismatch across the fixed cases is
6.14e-12 W. Particle-number changes remain below 7.4e-17 relative in the
independently reconstructed density-change integral. A tiny negative computed
vacuum energy (~1e-24 J) in the uniform control is constant-mode roundoff;
it is preserved, not clamped or presented as physical negative energy.

## Evolving-resistance response

The full seeded Spitzer run uses the same 8/16 grid and 320 steps. Final eta is
74.0080–74.0084 µΩ m and mean Te/Ti is 3.021323/2.978689 eV. Thermal gain is
2.6356634 µJ versus 2.4228428 µJ at fixed eta, an 8.7839% increase. Final maximum
current is 73.801 versus 350.083 A/m²; peak speed is 0.016232 m/s. These are
physical-model/control differences, not errors against the fixed-eta baseline.

Maximum energy error is 4.73137e-6 of the stated field scale. Peak electron Hall
parameter is 0.00599; current drift / electron and ion thermal speeds are
3.95e-5 and 0.01623. Those pass the existing review thresholds, not a proof of
complete closure validity. Conduction, thermal forces, viscosity, neutrals and
ionization remain absent from this conditional weak-field Z1 model.

## Evidence and next physical implementation

Actual sessions 66738 (uniform), 10198 (seeded coarse), 15544 (seeded fine),
13948 (time refinement) and 44808 (Spitzer) all exited 0 after their full runs.
Durations were 52.07, 53.29, 344.21, 90.77 and 167.67 s; peak RSS was 0.369,
0.375, 1.906, 0.518 and 0.389 GB. Plot observer 67929 also exited 0.
Outputs are `vacuum_pulse_fixed_r8_t160`, `vacuum_pulse_seed_r8_z16_t160`,
`vacuum_pulse_seed_r16_z32_t160`, `vacuum_pulse_seed_r8_z16_t320`, and
`vacuum_pulse_spitzer_r8_z16_t320`, under `outputs/`. Measured summaries and the
plot are in `outputs/vacuum_pulse_comparison`. The uniform case predates the
explicit nz option; both source versions are frozen. No worker remains live.

frc-vacuum.1 is complete only for this fixed-footprint periodic magnetic coupling.
NEXT frc-vacuum.2: connect finite circular-coil source fields through physical
axial end boundaries, preserving the current/force/heat and source-work coupling.
A periodic source array is not a substitute for those end effects. Exact moving
free surfaces, transport/ionization/validity, long-contact accuracy and converged
formation remain open. These runs establish neither laboratory recommendations,
FRC formation, three-dimensional stability nor particle confinement.
