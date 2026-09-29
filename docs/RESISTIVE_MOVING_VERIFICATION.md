# Resistive moving-plasma verification

The moving material/vacuum reference now includes scalar physical resistivity, consistent weak current, positive per-cell Ohmic heating and atomic Strang composition with ideal motion. Adaptive midpoint substeps control magnetic-field, current and local-heat error. No exterior fluid, pressure/density floor, numerical diffusivity or energy repair is introduced. The [numerical derivation](RESISTIVE_MOVING_NUMERICS.md) gives the equations, source-work contract and limits.

All 13 stationary diffusion cases and 20 moving cases are complete and source/configuration/field-hash checked. These are unit-scale mathematical verification inputs with rho=1 kg/m³ and gamma=5/3, without a chosen species/temperature/transport regime. The Stage 2 formation entry remains gated. Older ideal and fixed-chamber studies retain historical sources.

## Independent continuum diffusion, current and heat

A fixed conducting sphere of R=1 m has D=eta/mu0=0.1 m²/s and initial central a=A=sqrt(mu0) T. For k=n pi/R, a=3A j1(ks)/(ks) exp(-Dk²t) inside, with a continuously matched vacuum dipole outside. Jtheta=k² r a/mu0 inside and zero outside. Full-space energy is 6pi A²R/(mu0 k²) exp(-2Dk²t); integrated heat is its analytic loss. Independent Maxwell/Faraday/Ohm, interface and all-space quadrature tests validate this reference, separately from the numerical operator. See the [reference derivation](MOVING_RESISTIVE_PLAN.md).

The geometry is explicitly held fixed for this component test. Its magnetic force is not balanced by fluid pressure; this is not a freely evolving equilibrium. B is evaluated from actual Q1 gradients, recovered current from the consistent weighted mass solve, and local heat from integrated eta J². At t=0.1 s for the first radial mode:

| Radial cells | B volume L2 | Flux volume L2 | Current volume L2 | Local heat volume L2 | Boundary-current scaled max error | Axis-field relative error |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 0.09262226 | 0.03340527 | 0.05212094 | 0.124056 | 0.120361 | 0.03586362 |
| 8 | 0.04302237 | 0.008281384 | 0.01317189 | 0.02972155 | 0.02729424 | 0.01533338 |
| 16 | 0.02107254 | 0.002065239 | 0.003302659 | 0.007356166 | 0.006749845 | 0.005592423 |
| 32 | 0.0104809 | 0.0005160202 | 0.0008264589 | 0.00183494 | 0.001684015 | 0.001844669 |
| 64 | 0.005233511 | 0.000128996 | 0.0002066853 | 0.0004585813 | 0.0004208456 | 0.0005727495 |

The Q1-gradient B error is approximately first order; flux, current and local heat refine approximately at second order. Axis-field error has its own slower convergence. Boundary current is the material-side trace at the sphere vertices, normalized by the analytic equatorial current magnitude, not divided by a vanishing axis current. Local heat compares exact cell-integrated heat on actual polygonal cells. Global energy/heat also retain geometric error relative to the true sphere.

| Case | B relative L2 | Current relative L2 | Local heat relative L2 | Accepted diffusion substeps / rejected trials | Response residual [J] |
|---|---:|---:|---:|---:|---:|
| mesh_4 | 0.09262226 | 0.05212094 | 0.124056 | 100 / 0 | -6.217249e-15 |
| mesh_8 | 0.04302237 | 0.01317189 | 0.02972155 | 100 / 0 | 6.883383e-15 |
| mesh_16 | 0.02107254 | 0.003302659 | 0.007356166 | 100 / 0 | -1.554312e-15 |
| mesh_32 | 0.0104809 | 0.0008264589 | 0.00183494 | 102 / 1 | -2.220446e-15 |
| mesh_64 | 0.005233511 | 0.0002066853 | 0.0004585813 | 104 / 2 | -2.220446e-15 |
| time_half | 0.0104809 | 0.0008264641 | 0.001834895 | 200 / 0 | -3.330669e-15 |
| tolerance | 0.0104809 | 0.000826459 | 0.001834939 | 110 / 3 | 2.220446e-16 |
| far_matching | 0.02107999 | 0.003405325 | 0.007604285 | 100 / 0 | -1.332268e-15 |
| near_matching | 0.02105184 | 0.00300466 | 0.00664039 | 100 / 0 | -1.554312e-15 |
| modes4 | 0.02107254 | 0.003302324 | 0.007356269 | 100 / 0 | -3.108624e-15 |
| quadrature7 | 0.02107254 | 0.003302659 | 0.007356166 | 100 / 0 | -4.440892e-15 |
| mode2_16 | 0.0983412 | 0.005789365 | 0.009898087 | 104 / 2 | 2.775558e-16 |
| mode2_32 | 0.04924005 | 0.001448897 | 0.002493501 | 104 / 2 | -4.440892e-16 |

Time-half and tighter-tolerance controls use nr=32; matching, mode and quadrature controls use nr=16 with unchanged exterior radial spacing. `mode2` is the independently known second radial eigenmode, not numerically formed reversal. Its larger Q1-gradient field error remains visible even though current and heating errors are smaller. All configurations and entire sampled histories are archived.

## Stiff time integration: energy balance is insufficient

An independently diagonalized lowest spatial eigenmode gives an exact semidiscrete decay exp(-4). This isolates time integration from continuum spatial error. A single raw midpoint step returns amplitude -1/3: it falsely reverses the field while satisfying magnetic-loss/heat bookkeeping. The adaptive method retains positive decay and improves with tolerance:

| Local tolerance | Final amplitude / initial | Absolute amplitude error | Accepted substeps | Rejected trials |
|---:|---:|---:|---:|---:|
| 0.001 | 0.01798527863 | 0.0003303603 | 22 | 2 |
| 0.0001 | 0.01823359308 | 8.204581e-05 | 44 | 3 |
| 1e-05 | 0.01829663615 | 1.900274e-05 | 92 | 2 |
| 1e-06 | 0.01831151155 | 4.127343e-06 | 200 | 4 |

Exact amplitude is 0.01831563889; the raw-step balance residual is only 2.22e-16 J despite its incorrect sign. The controller accepts the uncorrected two-half-step solution and its physical heat. It does not add damping, clip the field or force an energy residual to zero. A local tolerance is not a global error or sign guarantee for arbitrarily small fields; topology claims need resolution and uncertainty controls.

## Actual moving resistive plasma

Each moving case begins in the analytic spherical FRC corresponding to its initial bias: a=-3B(0)(1-s²)/4 and p=(45/8)B(0)² r²(1-s²)/mu0. Positive eta then changes flux and deposits heat, so subsequent motion is resistive relaxation rather than pure equilibrium error. Source cases also apply B(t)=sqrt(mu0) exp(-(t-0.04)²/(2*0.04²)) T, ending at its peak at 0.04 s. Scalar D is prescribed mathematical physics, not a validated argon transport model.

| Case | D [m²/s] | Max speed [m/s] | Physical Ohmic heat [J] | Source work [J] | Response residual [J] | Outer steps | Diffusion substeps / rejects |
|---|---:|---:|---:|---:|---:|---:|---:|
| frc_4 | 0.1 | 0.07433695 | 0.3907056 | 0 | -7.254965e-09 | 80 | 320 / 0 |
| frc_8 | 0.1 | 0.07301204 | 0.3650503 | 0 | -2.832723e-08 | 80 | 320 / 0 |
| frc_16 | 0.1 | 0.1190891 | 0.3593711 | 0 | -3.134173e-08 | 140 | 560 / 0 |
| frc_32 | 0.1 | 0.1792189 | 0.3578613 | 0 | -7.043974e-09 | 531 | 2124 / 0 |
| source_4 | 0.1 | 0.1235151 | 0.3330767 | 1.9532285 | 0.0001114646 | 80 | 320 / 0 |
| source_8 | 0.1 | 0.1609598 | 0.3012837 | 1.9229666 | 0.0001097757 | 80 | 320 / 0 |
| source_16 | 0.1 | 0.1524294 | 0.287836 | 1.9103945 | 0.0001001414 | 90 | 360 / 0 |
| source_32 | 0.1 | 0.2657543 | 0.2832274 | 1.9066558 | 6.744483e-06 | 324 | 1296 / 0 |
| source_wide_16 | 1 | 0.3474634 | 1.591833 | 1.6081712 | 9.678204e-05 | 82 | 698 / 160 |
| source_wide_32 | 1 | 0.6669166 | 1.584777 | 1.609599 | 8.644224e-06 | 266 | 2018 / 416 |

Changing D from 0.1 to 1 m²/s changes physical diffusion/heating, and is not a numerical stabilization choice. The associated scale sqrt(D*0.04) is about 0.0632 or 0.2 m, respectively; this is only a diffusion-length estimate for a changing, moving geometry. The model has no exterior mass or Ohmic heat. Gas energy, field energy, physical heating and source work are separately recorded.

### Material-coordinate spatial comparisons

Corresponding material vertices supply velocity errors, using coarser inertial-mass weights. Actual Q1 fields/current are sampled at common reference-sector points with composite quadrature over every fine subcell. Density/pressure restrict by actual cell volume and integrated heat sums before division by volume. The two evolved coordinate maps can differ: these are material-coordinate state differences, not Cartesian pointwise errors. Changes subtract each mesh's initial state; total current is additionally compared without subtraction. All values below are maximum differences across the midpoint/final saved states, relative to the finer reference.

| Coarse / fine | Velocity [%] | Magnetic change [%] | Total current [%] | Current change [%] | Density change [%] | Pressure change [%] | Local heat density [%] |
|---|---:|---:|---:|---:|---:|---:|---:|
| frc_4 / frc_8 | 110.266 | 38.5724 | 8.22663 | 81.2378 | 126.84 | 19.1826 | 17.0037 |
| frc_8 / frc_16 | 49.2259 | 22.2043 | 2.09596 | 54.3286 | 40.8727 | 6.50299 | 4.44484 |
| frc_16 / frc_32 | 37.1209 | 11.8415 | 2.23379 | 26.97 | 35.2003 | 4.04015 | 1.78003 |
| source_4 / source_8 | 36.1607 | 55.666 | 22.9187 | 41.0597 | 38.0475 | 18.8154 | 20.3138 |
| source_8 / source_16 | 27.7989 | 49.5196 | 8.90699 | 16.1003 | 38.1971 | 10.064 | 11.1879 |
| source_16 / source_32 | 63.821 | 29.357 | 3.67623 | 8.02142 | 54.0889 | 6.07237 | 5.4924 |
| source_wide_16 / source_wide_32 | 53.9147 | 2.02999 | 2.12333 | 26.4342 | 48.6609 | 6.66961 | 4.35678 |

At D=0.1, the 16-to-32 velocity change is 63.821%; at D=1 it is 53.915%. These adjacent-grid differences are not absolute errors or proof of a fully converged formation experiment. Current/thermal/density changes must be judged alongside velocity and the fixed-sphere analytic checks.

### Time, tolerance and magnetic-boundary controls

The following use the same physical model and matching material resolution. A lower outer timestep does not replace diffusion error control, and a tighter diffusion tolerance does not replace material-wave/source-time resolution. Matching radius is varied at fixed exterior radial spacing.

| Test / reference | Velocity relative L2 [%] | Velocity RMS difference [m/s] | Magnetic-change L2 [%] | Total-current L2 [%] | Local-heat-density L2 [%] |
|---|---:|---:|---:|---:|---:|
| frc_time_0.001 / frc_8 | 0.004595021 | 1.181174e-06 | 0.001506729 | 0.0009199571 | 0.0003511726 |
| frc_8 / frc_time_0.00025 | 0.001170203 | 3.008087e-07 | 0.0003805267 | 0.0002469924 | 8.779855e-05 |
| frc_far_matching / frc_16 | 2.046628 | 0.0002937157 | 0.279068 | 0.03580376 | 0.06103979 |
| frc_tolerance / frc_16 | 0 | 0 | 0 | 0 | 0 |
| source_time_0.001 / source_8 | 0.05987686 | 1.048541e-05 | 0.006480152 | 0.004038692 | 0.02163062 |
| source_8 / source_time_0.00025 | 0.01496965 | 2.621165e-06 | 0.00161949 | 0.001010327 | 0.005407664 |
| source_far_matching / source_16 | 0.5397182 | 0.0002204233 | 0.07071837 | 0.0300722 | 0.0628732 |
| source_tolerance / source_16 | 0 | 0 | 0 | 0 | 0 |
| source_modes4 / source_16 | 0.002164077 | 8.02453e-07 | 0.001838516 | 0.0005646492 | 0.00134383 |
| source_quadrature7 / source_16 | 3.085865e-05 | 5.999778e-09 | 5.245204e-06 | 4.840387e-06 | 7.17832e-06 |

| Moving fixture | Outer time cap [s] | Steps | Response residual [J] | Ohmic heat [J] |
|---|---:|---:|---:|---:|
| frc, nr=8 | 0.001 | 40 | -1.1330805e-07 | 0.36505025 |
| frc, nr=8 | 0.0005 | 80 | -2.8327229e-08 | 0.36505028 |
| frc, nr=8 | 0.00025 | 160 | -7.0819048e-09 | 0.36505029 |
| source, nr=8 | 0.001 | 40 | 0.00043913151 | 0.3012901 |
| source, nr=8 | 0.0005 | 80 | 0.00010977569 | 0.30128372 |
| source, nr=8 | 0.00025 | 160 | 2.7443474e-05 | 0.30128213 |

## Interpretation and remaining requirements

All sampled material mass changes are 0 relative. Maximum sampled resistive-substep energy residual is 1.679e-14 J, while the full moving response has its separate finite-time error above. Maximum sampled vacuum relative residual is 4.451e-13. Minimum sampled density and pressure are 0.7366317 kg/m³ and 1.458604e-06 Pa. These sampled extrema are not bounds between outputs. Heat is nonnegative and actual cumulative cell heat is archived; there is no pressure floor.

Independent polynomial/volume checks validate the material-field sampler and conservative cell restrictions used in these comparisons. Atomic-failure tests cover a rejected second complete resistive half-step, after trial heating and material motion, and diffusion trial exhaustion. The accepted state, labels, entropy, heat, time and source work remain unchanged on failure.

This extends real moving-interface dynamics to physical resistive induction and local heating. It does not establish physical transport validity, surface thermodynamic closure, arbitrary coil geometry, shocks/remapping, long-time mesh robustness, Hall effects, ionization, particle confinement or FRC formation. The initial reversed-field topology is prescribed. Numerical field accuracy, current-layer resolution, actual boundary/thermal behavior and model validity remain independent requirements.

Next add spatially varying prescribed circular coils with the full geometric source-force derivative and source-work verification, while addressing the measured field/current/flow error and physical transport/Hall applicability. The dynamic experiment interface must preserve those validity limits rather than unlocking an unverified formation claim.

Full configurations, accepted histories, evolving labels/entropy, field/current/heat arrays, environment and hashes are in `outputs/resistive_moving_verification`. Reproduce with `python -m examples.resistive_decay_measurements`, `python -m examples.resistive_dynamic_measurements --reuse`, then `python -m examples.verify_resistive_moving`. The finest moving cases take minutes. Report/analysis/render sources are separately hashed from the frozen scientific measurements.

![Resistive accuracy](../outputs/resistive_moving_verification/resistive_verification.png)

![Computed resistive fields](../outputs/resistive_moving_verification/resistive_fields.png)
