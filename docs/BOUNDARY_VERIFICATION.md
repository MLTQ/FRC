# Boundary accuracy verification

This study tests a fixed first-order wall/end strip against the compatible second-order method. The shared conservative/CT assembly, physical resistivity, natural vacuum work and SSP RK3 are unchanged. The strip remains an explicit experimental probe; production defaults are unchanged. Fixed order reduction is recorded separately from adaptive fallback.

## Original shock, unchanged initial state

The independent central shock reference and the 0.02 s end time are identical to the compatible-energy study. Full-space response error remains distinct from exactly conservative fluid energy-flux accounting.

| Method | Radial cells | B L1/B1 | Pressure L1 [Pa] | Response residual [J] | Min density [kg/m³] | Additional adaptive limited stages |
|---|---:|---:|---:|---:|---:|---:|
| standard | 24 | 0.3094352 | 1.98346 | 2.971408e-05 | 0.28307 | 0 |
| standard | 48 | 0.2034171 | 1.121635 | 7.724324e-06 | 0.273282 | 0 |
| standard | 96 | 0.09993071 | 0.5912317 | 1.956309e-06 | 0.187014 | 24 |
| standard | 192 | 0.05931565 | 0.3343021 | 0.0005961242 | 0.0437523 | 550 |
| standard | 384 | 0.03065812 | 0.1753116 | 0.001732782 | 0.000916625 | 3423 |
| strip | 24 | 0.3095558 | 1.983687 | 2.467918e-05 | 0.285933 | 0 |
| strip | 48 | 0.2034179 | 1.121638 | 6.191567e-06 | 0.272084 | 0 |
| strip | 96 | 0.09993088 | 0.5912305 | 1.564342e-06 | 0.252714 | 0 |
| strip | 192 | 0.05931598 | 0.3343084 | 3.973689e-07 | 0.165396 | 0 |
| strip | 384 | 0.03065827 | 0.1753145 | 1.011297e-07 | 0.0496321 | 0 |

The nr=384 strip has 1534 forced cells and 3064 cells sharing selected faces/corners. Its adaptive count excludes these fixed selections. The earlier nr=192 time control decreases response residual from 3.973689e-7 to 4.788792e-8 J on CFL halving. The new fine-grid result alone does not replace a separate nr=384 temporal study.

## Independent stationary-cylinder skin test

For an infinite cylinder, Bz solves B_t=D(B_rr+B_r/r), with prescribed exterior Bz=Re(Bs exp(iωt)), regular axis and continuous boundary B. The independent solution uses k=(1+i)/δ, δ=√(2D/ω), Bhat=Bs I0(kr)/I0(kR), psihat=Bs r I1(kr)/(k I0(kR)), and Jhat=-Bs k I1(kr)/(μ0 I0(kR)). Outside, J=0 and psi continues as psi(R)+Bs(r²-R²)/2. The exact periodic solution initializes every case; eight phases over one period are measured.

The magnetic-diffusion and skin-effect basis follows [MIT, magnetic diffusion laws](https://web.mit.edu/6.013_book/www/chapter10/10.5.html) and [MIT, skin effect](https://www.mit.edu/course/6/6.013_book/www/chapter10/10.7.html). The cylindrical solution is derived here by separation. Independent finite differences check Maxwell/Ohm/diffusion identities, and Gauss volume integration of mean ηJ² agrees with analytic inward surface Poynting power.

This is a kinematic conductor check: gas is held at rho=1 kg/m³, p=1 Pa and v=0 to supply the declared production LLF speed. Standard/strip cases use the actual production reconstruction, CT electric fields and wall numerical resistance. Only the invariant-z central stencil is evolved; end/corner dynamics and gas motion are not tested here. The physical-only control omits both UCT and numerical interface diffusion. It is not a full MHD option.

| Method | Radial cells | δ / grid spacing | Maximum B L2 / Bs | Maximum wall J error / (Bs/μ0δ) | Maximum flux error / (Bs R²/2) |
|---|---:|---:|---:|---:|---:|
| physical | 16 | 0.9027 | 0.06875587 | 0.3642259 | 0.02619825 |
| physical | 32 | 1.8054 | 0.01952049 | 0.1071207 | 0.007330655 |
| physical | 64 | 3.6108 | 0.004857269 | 0.0270297 | 0.001841247 |
| physical | 128 | 7.2216 | 0.001209741 | 0.006756788 | 0.0004599448 |
| standard | 16 | 0.9027 | 0.05071535 | 0.3321142 | 0.01993159 |
| standard | 32 | 1.8054 | 0.01426935 | 0.1226187 | 0.004963483 |
| standard | 64 | 3.6108 | 0.00367116 | 0.05904579 | 0.001271504 |
| standard | 128 | 7.2216 | 0.0009344137 | 0.03203018 | 0.0003209786 |
| strip | 16 | 0.9027 | 0.05335723 | 0.3731367 | 0.02513184 |
| strip | 32 | 1.8054 | 0.01727821 | 0.09823902 | 0.007335766 |
| strip | 64 | 3.6108 | 0.005186479 | 0.03247238 | 0.002306257 |
| strip | 128 | 7.2216 | 0.001507676 | 0.02407089 | 0.0007036572 |

R=1 m, D=0.1 m²/s, f=10 Hz give δ=0.05641896 m. B errors compare exact cylindrical cell averages; norms use fixed applied-field scales, not a near-zero instantaneous norm. Mean analytic physical dissipation is 5.41014368e-08 W/m for the recorded tiny field amplitude. Surface current has a first-order pointwise stencil and is separately measured. Partial cancellation can make a numerically diffusive result closer than physical-only on a particular grid; that is not additional physical transport.

| Method, nr=64 | CFL 0.2 maximum B error | CFL 0.1 maximum B error | Relative change [%] |
|---|---:|---:|---:|
| physical | 0.0048572686 | 0.0048572874 | -0.00038805 |
| standard | 0.0036711595 | 0.0036711549 | 0.00012656 |
| strip | 0.0051864785 | 0.0051864591 | 0.00037489 |

## Identical finite-cylinder coil pulses

Both methods use the original two-coil pulse, fixed gas chamber, infinite multipole magnetic exterior, physical D=0.1 m²/s and 0.04 s rise. Source amplitude, density, pressure, and all constitutive/numerical controls are recorded. Full conserved, primitive, CT and current fields are saved at 0.04, 0.06 and 0.10 s in lossless hashed NPZ files.

| Case | Source work [J] | Response residual / work [%] | Final peak speed [m/s] | Final kinetic energy [J] | Numerical / physical heat [%] | Additional adaptive limited stages |
|---|---:|---:|---:|---:|---:|---:|
| standard_16 | 0.130062259 | 0.02992537 | 0.129674 | 0.007396138 | 18.0189 | 0 |
| standard_32 | 0.125537166 | 0.01096728 | 0.179401 | 0.01100596 | 6.00842 | 0 |
| standard_64 | 0.121832577 | 0.00315436 | 0.2075944 | 0.01160384 | 1.78759 | 0 |
| standard_128 | 0.120476 | 0.0008375808 | 0.2218445 | 0.0118463 | 0.480364 | 0 |
| strip_16 | 0.143408696 | 0.03044819 | 0.1265123 | 0.00656207 | 20.6058 | 0 |
| strip_32 | 0.130398402 | 0.01060995 | 0.1673018 | 0.008685492 | 6.66667 | 0 |
| strip_64 | 0.123820317 | 0.003113323 | 0.1997911 | 0.01054832 | 1.86701 | 0 |
| strip_128 | 0.121272192 | 0.0008346501 | 0.2178213 | 0.01150044 | 0.490463 | 0 |
| standard_half_time | 0.121832576 | 0.003153454 | 0.2075944 | 0.01160384 | 1.78759 | 0 |
| strip_half_time | 0.123820315 | 0.003112394 | 0.1997911 | 0.01054832 | 1.86701 | 0 |
| standard_far_matching | 0.125537812 | 0.01096723 | 0.1794014 | 0.01100603 | 6.00842 | 0 |
| strip_far_matching | 0.130399069 | 0.01060989 | 0.1673022 | 0.008685543 | 6.66668 | 0 |

Source work excludes external-field self/mutual and coil-circuit energy. Numerical interface heat is separate from physical Ohmic heat and does not include all implicit Rusanov/UCT dissipation. The characteristic rise diffusion length √(D t_rise)=0.0632456 m spans about 1.01, 2.02, 4.05 and 8.10 radial cells on these grids; this is a scale estimate, not the exact finite-coil skin profile.

| Method at nr=128 | Time [s] | Density min / max [kg/m³] | J²-weighted RMS depth [m] | Depth / radial cell | Max local numerical / physical eta [%] |
|---|---:|---:|---:|---:|---:|
| standard | 0.04 | 0.807673 / 1.05094 | 0.0441975 | 5.6573 | 4.31867 |
| standard | 0.06 | 0.72497 / 1.09281 | 0.0703355 | 9.0029 | 4.40427 |
| standard | 0.1 | 0.816407 / 1.12608 | 0.102789 | 13.157 | 3.48023 |
| strip | 0.04 | 0.868175 / 1.04549 | 0.0447746 | 5.7311 | 4.28758 |
| strip | 0.06 | 0.785843 / 1.08781 | 0.0714223 | 9.1421 | 4.32171 |
| strip | 0.1 | 0.812954 / 1.12243 | 0.102911 | 13.173 | 3.43349 |

Current depth is measured from the nearest gas wall/end with cylindrical nodal J² weights. It is a width diagnostic rather than an assumed exponential fit. Local eta ratios and integrated heat ratios measure different aspects of numerical dissipation.

### Full-field spatial differences

Cell fields are restricted from the finer mesh by exact cylindrical volume weights; currents use coincident nodes with cylindrical dual weights. Values below are the largest relative L2 difference across the three saved times. Density/pressure/B denominators use perturbations from the initial gas/bias; velocity/current use their reference norms. These are adjacent-grid differences, not analytic errors or proof of the finest grid's absolute accuracy.

| Coarse / fine | Velocity [%] | Density perturbation [%] | Pressure perturbation [%] | Magnetic perturbation [%] | Current [%] |
|---|---:|---:|---:|---:|---:|
| standard_16 / standard_32 | 21.4364 | 35.5737 | 53.0438 | 9.75575 | 10.7932 |
| standard_32 / standard_64 | 21.8188 | 30.7604 | 41.31 | 5.67117 | 11.2346 |
| standard_64 / standard_128 | 7.05911 | 12.1023 | 12.8073 | 1.04401 | 5.36545 |
| strip_16 / strip_32 | 23.8263 | 46.9162 | 48.2248 | 12.6525 | 16.2574 |
| strip_32 / strip_64 | 18.2401 | 44.9229 | 42.511 | 5.31381 | 12.9801 |
| strip_64 / strip_128 | 10.4227 | 34.1558 | 33.5496 | 2.47763 | 9.3832 |

### Same-grid method, time and matching controls

The time control halves both CFL and the prescribed-source time cap. The farther matching control preserves the same actual gas domain and sources.

| Test / reference | Max velocity L2 difference [%] | Max magnetic perturbation L2 difference [%] | Max current L2 difference [%] |
|---|---:|---:|---:|
| strip_16 / standard_16 | 28.7522 | 11.1112 | 22.657 |
| strip_32 / standard_32 | 24.5915 | 5.81165 | 12.674 |
| strip_64 / standard_64 | 15.9224 | 2.61017 | 10.2022 |
| strip_128 / standard_128 | 9.76713 | 1.18472 | 5.75412 |
| standard_64 / standard_half_time | 7.09866e-05 | 2.42515e-05 | 0.000246754 |
| standard_32 / standard_far_matching | 0.000332145 | 0.000271724 | 0.000331886 |
| standard_64 / standard_length1_64 | 7.42592e-05 | 6.02172e-05 | 7.43156e-05 |
| strip_64 / strip_half_time | 2.34586e-05 | 2.44536e-05 | 0.000252987 |
| strip_32 / strip_far_matching | 0.000321613 | 0.000270975 | 0.000327036 |
| strip_64 / strip_length1_64 | 7.30364e-05 | 6.00519e-05 | 7.38769e-05 |

## Gas-end sensitivity with a common magnetic boundary

At nr=64 and matching extent 3 m, move reflecting gas ends from z=±1 to ±2 m with identical coils, initial gas state and spacing. The physical gas domain and total mass change. Compare common cells at |z|<0.5 m; source work refers to the whole actual domain. This quantifies a boundary-condition choice rather than demonstrating an open exhaust or free surface.

| Method | Short/long source-work difference [% of long] | Max central velocity difference [%] | Max central magnetic perturbation difference [%] |
|---|---:|---:|---:|
| standard | 2.10556 | 0.817383 | 0.960513 |
| strip | 2.118 | 0.837169 | 0.977936 |

## Interpretation and remaining scope

The strip has a material accuracy cost in the driven gas: at nr=128 the maximum same-grid velocity L2 difference is 9.767%. The 64-to-128 changes are 7.059% for the standard method and 10.42% for the strip. These measured differences prevent treating the strip's small shock energy residual as evidence of an accurate general-purpose boundary replacement. It stays experimental and no production default is promoted.

The fixed strip is a controlled numerical option, not material physics. Assess its shock energy and analytic skin results together with driven-flow grid differences and the changed-end controls. A positive state, small response residual, or close agreement between two methods cannot establish a physically accurate formation experiment.

The evolving free-plasma/vacuum interface, actual gas loss/ends, regime-appropriate transport and ionization, broader shock/cold-flow and curved-equilibrium accuracy, converged formation and dynamic experiment UI remain open. No exterior fluid atmosphere, energy correction or transport closure is introduced by this study. Stage 2 formation remains gated.

Full inputs, independent references, profiles, field hashes and source/environment records are in `outputs/boundary_verification/metrics.json` and its case directories. The original compatible-energy records retain their own provenance. Rebuild with `python -m examples.verify_boundary` after all cases finish.

![Boundary accuracy study](../outputs/boundary_verification/boundary_verification.png)
