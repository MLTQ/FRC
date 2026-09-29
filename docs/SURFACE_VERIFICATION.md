# Resolving the moving resistive surface layer

This study adds explicit radial grading to the actual material/vacuum mesh while retaining the same physical D=eta/mu0=0.1 m²/s. Material shells use R[1-(1-xi)^p], and vacuum offsets use (Rmatch-R)xi^p. The default p=1 retains uniform spacing. This is an initial mesh choice, not adaptive remapping, artificial resistance or a transport closure.

All 21 moving cases and 13 independent component cases are complete with source/configuration/field hashes checked. The extra finest-mesh time control was added after the first control exposed time sensitivity, using the same frozen measurement entry point and fully recorded typed inputs; its reproduction helper is separately hashed. The cases use the preceding mathematical FRC and uniform-bias fixtures, rho=1 kg/m³ and gamma=5/3. Reversed flux is initialized. These are not demonstrated formation experiments or validated species parameters.

## Independent graded-mesh checks

Physical volume and second moments match independent disk integration of the revolved polygon. Actual free-gas expansion is checked against lambda=sqrt(1+2t²), with p0=1-s² and rho0=1. Fixed conducting-sphere decay independently checks full profiles, current, local Ohmic heat and surface current. The gas test has no heat source; the diffusion test has fixed geometry. Neither alone validates a heated free surface.

| Decay case | B relative L2 | Current relative L2 | Local heat relative L2 | Surface-current scaled maximum error |
|---|---:|---:|---:|---:|
| decay_g1.5_8 | 0.04246511 | 0.01130405 | 0.0250676 | 0.02301439 |
| decay_g1.5_16 | 0.02086087 | 0.002840129 | 0.006258987 | 0.005724227 |
| decay_g1.5_32 | 0.01038371 | 0.0007113753 | 0.001565729 | 0.001430254 |
| decay_g2_8 | 0.04298661 | 0.01313455 | 0.02903657 | 0.02503133 |
| decay_g2_16 | 0.02110838 | 0.003308105 | 0.007270903 | 0.006229101 |
| decay_g2_32 | 0.01050699 | 0.0008292687 | 0.001820607 | 0.001556561 |
| decay_g2_16_angular | 0.0120315 | 0.00284562 | 0.006010668 | 0.003980721 |

| Free-gas case, t=0.04 s | Velocity relative mass L2 | Maximum surface-position error [m] | Response-energy residual [J] |
|---|---:|---:|---:|
| expansion_g1.5_4 | 0.1128126 | 0.0003291985 | -3.409499e-09 |
| expansion_g1.5_8 | 0.04766092 | 0.0002327912 | -3.134357e-09 |
| expansion_g1.5_16 | 0.01923881 | 0.000160808 | -3.041682e-09 |
| expansion_g2_4 | 0.1589823 | 0.0006439818 | -4.150364e-09 |
| expansion_g2_8 | 0.07561347 | 0.0006307053 | -3.791555e-09 |
| expansion_g2_16 | 0.03075555 | 0.0005665697 | -3.356555e-09 |

Current and local heat retain approximately second-order refinement, while Q1-gradient B remains about first order. Grading trades interior resolution for surface resolution; the smooth gas-expansion error therefore also needs measurement. No automatic accuracy claim follows from a smaller last cell.

## Actual moving fields and surface motion

All moving runs end at 0.04 s. Source cases use B(t)=sqrt(mu0) exp(-(t-0.04)²/(2*0.04²)) T; unforced FRC cases hold B=sqrt(mu0) T. Positive resistivity changes flux and deposits heat. The initial pressure vanishes at the analytic surface while current does not, so surface thermal evolution is a distinct physical and numerical requirement.

| Case | nr × nt | Grading | Peak speed [m/s] | Equatorial outward surface speed [m/s] | Minimum density [kg/m³] | Physical heat [J] | Steps |
|---|---:|---:|---:|---:|---:|---:|---:|
| source_g1_8 | 8 × 16 | 1 | 0.1609598 | -0.1609598 | 0.9929567 | 0.3012837 | 80 |
| source_g1_16 | 16 × 32 | 1 | 0.1524294 | -0.1518424 | 0.9925789 | 0.287836 | 90 |
| source_g1_32 | 32 × 64 | 1 | 0.2657543 | -0.04495779 | 0.9479444 | 0.2832274 | 324 |
| source_g1.5_8 | 8 × 16 | 1.5 | 0.1858897 | -0.08134645 | 0.9940058 | 0.2896935 | 80 |
| source_g1.5_16 | 16 × 32 | 1.5 | 0.3063933 | 0.1685147 | 0.7140144 | 0.2826868 | 80 |
| source_g1.5_32 | 32 × 64 | 1.5 | 0.3749719 | 0.3458967 | 0.3055524 | 0.2786056 | 220 |
| source_g2_8 | 8 × 16 | 2 | 0.2891885 | 0.2341573 | 0.6802363 | 0.2880219 | 80 |
| source_g2_16 | 16 × 32 | 2 | 0.4269923 | 0.3815274 | 0.2236513 | 0.2801046 | 84 |
| source_g2_32 | 32 × 64 | 2 | 0.5834666 | 0.1310004 | 0.0565282 | 0.2750341 | 210 |
| frc_g1.5_16 | 16 × 32 | 1.5 | 0.2880509 | 0.2880509 | 0.7737698 | 0.3583303 | 100 |
| frc_g1.5_32 | 32 × 64 | 1.5 | 0.3806321 | 0.3806321 | 0.4584842 | 0.3565137 | 358 |
| frc_g2_16 | 16 × 32 | 2 | 0.408497 | 0.408497 | 0.3401889 | 0.3575176 | 98 |
| frc_g2_32 | 32 × 64 | 2 | 0.4092852 | 0.3726458 | 0.09754458 | 0.3551418 | 307 |
| source_g2_16_angular | 16 × 64 | 2 | 0.4262247 | 0.3802864 | 0.2250495 | 0.2800327 | 92 |
| source_g2_32_angular | 32 × 128 | 2 | 0.5780913 | 0.1110813 | 0.05609764 | 0.274364 | 327 |
| source_g2_16_time | 16 × 32 | 2 | 0.4341877 | 0.3685921 | 0.2128477 | 0.2791347 | 160 |
| source_g2_32_cfl | 32 × 64 | 2 | 0.5727181 | 0.1169288 | 0.0557784 | 0.2741547 | 406 |
| source_g2_16_tolerance | 16 × 32 | 2 | 0.4269738 | 0.3815085 | 0.2236577 | 0.2801052 | 84 |
| source_g2_16_matching | 16 × 32 | 2 | 0.4180131 | 0.3731748 | 0.2279215 | 0.2795845 | 84 |
| source_g2_16_quadrature | 16 × 32 | 2 | 0.4269876 | 0.381509 | 0.2236456 | 0.2801042 | 84 |
| source_g2_32_cfl_quarter | 32 × 64 | 2 | 0.5668651 | 0.1214068 | 0.055544 | 0.2738125 | 798 |

The surface and the immediately interior material can move in opposite directions. The current, density, pressure and complete surface profiles are archived and plotted; bulk flow norms must not conceal a poorly resolved boundary. No floor, exterior gas or energy repair is introduced.

## Comparison on common initial physical material coordinates

The sampler inverts each initial straight-sector map, including its graded radial shells, and evaluates the actual Q1 fields at the same initial (r,z). Velocity uses initial physical mass weights; magnetic/current/thermal comparisons use evolved test-volume weights. Field, density and pressure changes subtract each mesh's own initial sample. This differs from the prior coarse-vertex velocity norm, so the percentages below must not replace or be equated to that report's 63.821% number. Thermal fields also differ in definition: the earlier report conservatively restricted fine-cell values to coarse cells; this pointwise P0 comparison retains differences inside those coarse cells.

The angular-coarse polygon lies inside the angular-fine polygon. Its support omits a reported fraction of the finer initial volume. A separate maximum boundary-vertex velocity difference compares shared initial surface locations and retains the surface-motion uncertainty. Differences from a finer run are not absolute continuum errors.

| Test / reference | Velocity mass L2 [%] | Magnetic change [%] | Total current [%] | Pressure change [%] | Density change [%] | Local heat density [%] | Max boundary-speed difference [m/s] | Omitted initial volume [%] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| source_g1_8 / source_g1_16 | 23.3033 | 50.488 | 10.4038 | 46.5224 | 48.6008 | 47.4668 | 0.01604597 | 0.72171 |
| source_g1_16 / source_g1_32 | 34.1273 | 29.9195 | 3.75219 | 31.9964 | 94.5936 | 33.9855 | 0.1068846 | 0.180645 |
| source_g1.5_8 / source_g1.5_16 | 47.8581 | 27.1934 | 10.0597 | 28.2754 | 99.2421 | 30.7358 | 0.2498612 | 0.72171 |
| source_g1.5_16 / source_g1.5_32 | 46.4582 | 13.8734 | 12.2861 | 17.4966 | 92.5633 | 15.6031 | 0.2301489 | 0.180645 |
| source_g2_8 / source_g2_16 | 63.3015 | 23.7857 | 17.7867 | 23.174 | 121.168 | 24.839 | 0.2474297 | 0.72171 |
| source_g2_16 / source_g2_32 | 55.2917 | 14.2694 | 17.8201 | 24.8438 | 148.338 | 13.0976 | 0.252004 | 0.180645 |
| frc_g1.5_16 / frc_g1.5_32 | 50.5653 | 10.096 | 4.52575 | 7.30465 | 81.628 | 5.24074 | 0.1027255 | 0.180645 |
| frc_g2_16 / frc_g2_32 | 77.7765 | 10.7076 | 7.08262 | 8.19584 | 105.975 | 4.54742 | 0.1122273 | 0.180645 |
| source_g1_32 / source_g1.5_32 | 61.2152 | 17.7064 | 14.4096 | 20.9894 | 100.424 | 21.8208 | 0.4042947 | 0 |
| source_g1.5_32 / source_g2_32 | 53.4682 | 13.1745 | 15.9755 | 23.3725 | 111.424 | 12.5516 | 0.310869 | 0 |
| source_g2_16 / source_g2_16_angular | 22.6825 | 6.2849 | 4.74578 | 10.1458 | 45.7375 | 6.22651 | 0.01234741 | 0.180645 |
| source_g2_16 / source_g2_16_time | 4.89221 | 0.883274 | 1.3543 | 1.72876 | 6.30685 | 1.49351 | 0.02626768 | 0 |
| source_g2_16 / source_g2_16_tolerance | 0.0022724 | 0.000759452 | 0.00542533 | 0.000970154 | 0.00279737 | 0.00107362 | 2.070021e-05 | 0 |
| source_g2_16 / source_g2_16_matching | 1.77707 | 0.262977 | 0.313573 | 0.182914 | 2.58447 | 0.166729 | 0.009230137 | 0 |
| source_g2_16 / source_g2_16_quadrature | 0.00356034 | 0.00182759 | 0.0143728 | 0.00105655 | 0.00330251 | 0.00087251 | 1.968685e-05 | 8.88178e-14 |
| source_g2_32 / source_g2_32_angular | 15.8362 | 2.93997 | 4.23203 | 11.8019 | 34.8629 | 6.43548 | 0.03575655 | 0.0451749 |
| source_g2_32 / source_g2_32_cfl | 5.39874 | 1.37631 | 1.17915 | 2.02683 | 4.17165 | 1.81537 | 0.04457785 | 0 |
| source_g2_32_cfl / source_g2_32_cfl_quarter | 2.758 | 0.721877 | 1.26636 | 0.997303 | 2.60914 | 0.881899 | 0.0181078 | 0 |

Each table entry is the maximum across the midpoint/final saved times. Pressure/density/heat remain actual discontinuous cell fields; no smoothing is used. Initial material-cell intersections form a common integration mesh. Positive triangle rules integrate inside each overlap, and actual deformed Jacobian/radius ratios provide volume weights. The integrated overlap volume is checked against the whole test polygon. Comparison quadrature is checked independently by increasing its one-dimensional Gauss order from 4 to 6:

| Norm | Maximum relative change in the measured difference [%] |
|---|---:|
| velocity_relative_mass_L2 | 3.27207e-10 |
| total_current_relative_volume_L2 | 1.12027e-08 |
| local_heat_density_relative_volume_L2 | 6.46935e-14 |
| magnetic_change_relative_volume_L2 | 8.84888e-09 |
| current_change_relative_volume_L2 | 3.37344e-07 |
| pressure_change_relative_volume_L2 | 2.82972e-14 |
| density_change_relative_volume_L2 | 9.55362e-14 |

These integration sensitivities are part of the diagnostic uncertainty, especially for piecewise constant thermal fields. They are not time-step errors. Explicit time/CFL, angular, diffusion tolerance, magnetic-matching and physical quadrature controls are distinct rows above.

The finest grade-2 mesh still changes by 15.8362% in velocity when angular resolution doubles. Halving CFL and the outer step cap gives 5.39874%, then 2.75800% velocity differences on the same mesh. Their ratio is about 1.96, not evidence of second-order asymptotic convergence; current differences do not decrease monotonically. These controls retain their fixed diffusion tolerance. Surface-speed differences remain 0.04457785 then 0.01810780 m/s. No time-extrapolated solution or resolved boundary is claimed.

## Force cancellation and accounting

The following final equatorial surface accelerations come from restoring the saved accepted state. The recovered-volume J cross B column separately projects weak current times Q1 B. Both magnetic columns are discrete constructions; their agreement is a diagnostic, not an independent continuum-force proof.

| Case | Pressure [m/s²] | Magnetic virtual work [m/s²] | Recovered volume J cross B [m/s²] | Net [m/s²] |
|---|---:|---:|---:|---:|
| source_g1_16 | 11.62341 | -12.53834 | -11.97937 | -0.9149315 |
| source_g1.5_32 | 53.18911 | -54.80362 | -53.94303 | -1.61451 |
| source_g2_16 | 72.67753 | -74.01166 | -74.58659 | -1.334131 |
| source_g2_32 | 259.4958 | -264.4388 | -261.8353 | -4.943024 |
| source_g2_32_angular | 264.0266 | -264.3109 | -261.7103 | -0.2843352 |

Maximum sampled relative mass change is 0; resistive-stage energy residual is 1.386e-13 J; vacuum relative residual is 7.566e-13. Whole moving-response time error is independently retained in every history and is not repaired. Physical heat is nonnegative in all saved cells. These sampled statements are not bounds between outputs.

The new uniform nr=8/16/32 runs have bitwise identical field/state arrays compared with the prior verified uniform-grid checkpoint. Historical results keep their historical fingerprints.

## Remaining scope

A separate initial-time compatibility audit derives a positive equatorial Joule pressure source of 3.75 Pa/s (unforced FRC) and 1.379547904 Pa/s (initial source fixture). If B is continuous, no additional surface stress/heat term is present and p=0 remains a smooth material trace with bounded velocity gradient, its material time derivative would instead be zero. Thus these initial heated-edge data fail that conditional smooth-boundary compatibility requirement. Independent one-sided Maxwell derivatives check the current used in the calculation. See [the derivation, assumptions and primary sources](SURFACE_COMPATIBILITY.md). This does not prove nonexistence of every weak solution or excuse a numerical error; it changes how these initialized resistive runs can be used as convergence targets.

Surface resolution is now an explicit, independently checked experiment control. The measured changes and boundary-speed differences remain the evidence of what is and is not resolved. This does not establish heated-surface continuum closure, shock/remapping robustness, temperature-dependent transport, ionization or Hall validity, actual chamber/end boundaries, particle confinement or FRC formation. Spatial coil coupling remains the next source capability, with its independent force reference in MOVING_COIL_PLAN.md. Stage 2 remains gated.

Reproduce the component records with `python -m examples.surface_measurements --references`, the original 20 moving records with `python -m examples.surface_measurements --reuse`, the additional time control with `python -m examples.surface_time_control`, and this assembly with `python -m examples.verify_surface`. The finest angular/time controls take many minutes. All source/configuration/field hashes and separate analysis/render/reproduction sources are retained in `outputs/surface_verification/metrics.json`.

![Computed surface profiles](../outputs/surface_verification/surface_profiles.png)

![Surface accuracy](../outputs/surface_verification/surface_accuracy.png)
