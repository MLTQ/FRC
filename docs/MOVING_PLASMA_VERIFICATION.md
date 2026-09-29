# Moving material plasma / infinite vacuum verification

This study verifies an experimental ideal, adiabatic material mesh with an actual moving plasma surface and no exterior fluid. A weighted regular-axis magnetic finite-element solve is matched continuously to analytic decaying modes at infinity. The [method derivation](MOVING_PLASMA.md) states equations, quadrature, source work, advection approximations and limitations. Existing fixed-chamber methods and defaults are unchanged; Stage 2 formation remains gated.

The 26 evolution cases and 27 static evaluations use SI unit-scale mathematical fixtures, not a selected gas/temperature/transport regime. Every evolution begins independently. All complete configurations, accepted histories, source/environment hashes and lossless field archives are retained in `outputs/moving_verification`. Old fixed-boundary archives retain their own historical sources.

## Independent magnetic energy and generalized force

For a unit initial sphere with Bp0=sqrt(mu0), uniform bias=0.3 Bp0 and homologous material scale lambda, fixed flux gives Bp=Bp0/lambda². Integrating the interior uniform induced field and the exterior dipole gives E=pi(1/lambda-0.6 lambda+0.09 lambda³) J. Its independent derivative is pi(-1/lambda²-0.6+0.27 lambda²). No assembled operator is used to construct these targets. This electromagnetic fixture is not a valid free-gas equilibrium with its nonzero normal field and tangential-field jump.

| Scale | Radial cells | Energy relative error | Generalized-force relative error |
|---|---:|---:|---:|
| 0.9 | 4 | 0.0425438 | 0.1547178 |
| 0.9 | 8 | 0.01132491 | 0.04100617 |
| 0.9 | 16 | 0.002882172 | 0.01042201 |
| 0.9 | 32 | 0.0007240942 | 0.00261671 |
| 1.0 | 4 | 0.01897219 | 0.09027764 |
| 1.0 | 8 | 0.005115946 | 0.02381069 |
| 1.0 | 16 | 0.001307549 | 0.0060408 |
| 1.0 | 32 | 0.0003290268 | 0.001516033 |
| 1.1 | 4 | 0.002781274 | 0.04784543 |
| 1.1 | 8 | 0.0008851846 | 0.01262142 |
| 1.1 | 16 | 0.000237522 | 0.003202816 |
| 1.1 | 32 | 6.07605e-05 | 0.0008040054 |

The finite material sphere is polygonal and has second-order geometric error; the distant matching sphere is exactly curved. Force is the analytic derivative of the full discrete potential, including vacuum mesh geometry, material flux coefficients and the regular-axis chain rule. An arbitrary two-dimensional deformation is separately checked against a symmetric perturbation of that complete potential. Agreement of two discrete energy quadratures does not remove continuum error.

### Translation, modes and matching

Axial translation in uniform bias leaves exact full-space energy unchanged. The finite magnetic mesh breaks this symmetry slightly; its net axial force is a numerical error, not a physical centering force.

| Radial cells | Shift [m] | Spurious axial force [N] | Energy relative error |
|---:|---:|---:|---:|
| 4 | 0.3 | -0.007965873 | 0.01972608 |
| 8 | 0.3 | -0.001970671 | 0.005302486 |
| 16 | 0.3 | -0.0004909631 | 0.001354038 |
| 32 | 0.3 | -0.0001223934 | 0.0003406183 |

The next controls use 12 radial cells and a 0.4 m axial translation. Moving the matching sphere preserves exterior radial spacing, but changes exterior element shapes; matching sensitivity can persist even with a tiny elliptic residual.

| Control | Value | Energy relative error | Generalized-force relative error | Axial force [N] |
|---|---:|---:|---:|---:|
| modes | 1 | 0.01983308 | 0.007581186 | -0.1503495 |
| modes | 2 | 0.00381969 | 0.01066714 | -0.02274863 |
| modes | 4 | 0.002467987 | 0.01084687 | -0.00141282 |
| modes | 8 | 0.002461094 | 0.01084698 | -0.001222816 |
| modes | 12 | 0.002461064 | 0.01084698 | -0.001222508 |
| matching | 1.5 | 0.002152842 | 0.01391367 | -0.00689161 |
| matching | 2 | 0.002461094 | 0.01084698 | -0.001222816 |
| matching | 3 | 0.002632182 | 0.008811318 | -0.000283103 |
| quadrature | 3 | 0.002461023 | 0.01084706 | -0.001222583 |
| quadrature | 5 | 0.002461094 | 0.01084698 | -0.001222816 |
| quadrature | 7 | 0.002461094 | 0.01084698 | -0.001222816 |

Maximum stationary-vacuum relative residual in the static study is 1.241e-13. Increasing modes eventually reaches spatial-error limits. Exact magnetic axial-momentum conservation is not claimed.

## Actual free gas expansion

For gamma=5/3, rho0=1 kg/m³, p0=1-s² Pa, B=0 and zero initial velocity inside the unit sphere, lambda''=2/lambda³ gives lambda=sqrt(1+2t²), x=lambda X, v=2t X/lambda. Surface pressure is zero and the exterior has no mass. The table compares actual evolution to this independent solution at 0.1 s. Velocity uses positive initial nodal-mass weights; pressure compares exact cell averages on the actual final geometry.

| Radial cells | Steps | Velocity relative L2 | Maximum surface-position error [m] | Density relative L2 | Pressure cell-average relative L2 | Energy residual [J] |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 100 | 0.07263601 | 0.001059381 | 0.00406602 | 0.003278419 | -7.600901e-08 |
| 8 | 100 | 0.02682206 | 0.0005611963 | 0.002684422 | 0.00112576 | -7.259995e-08 |
| 16 | 300 | 0.009249362 | 0.0002833768 | 0.00175235 | 0.0003571342 | -8.301605e-09 |
| 32 | 1147 | 0.003022045 | 0.0001382236 | 0.0011269 | 0.0001016082 | -5.487362e-10 |

The exact final radius is 1.009950494 m. Initial polygonal volume error and P0 entropy projection are retained; no correction makes the reference exact on a coarse mesh. Refinement of the velocity and surface errors provides stronger evidence than exact material mass alone.

## Initialized spherical FRC evolution

Initialize the independently derived [spherical FRC](FREE_BOUNDARY_REFERENCE.md) with B0=sqrt(mu0), unit radius, rho=1 kg/m³ and gamma=5/3. Inside, a=-3 B0(1-s²)/4 and p=(45/8)r²(1-s²) Pa. Outside, uniform bias plus the solved decaying field supplies vacuum traction. The continuum state has zero acceleration; all computed motion below is numerical drift. Its exact internal and induced energies are 9pi/7 and 10pi/7 J.

| Radial cells | Initial mass-RMS acceleration [m/s²] | Max speed at 0.04 s [m/s] | Max surface displacement [m] | Energy residual [J] |
|---:|---:|---:|---:|---:|
| 4 | 0.7130865 | 0.07642413 | 0.0006279655 | -2.217166e-08 |
| 8 | 0.2433248 | 0.03278417 | 0.000421427 | -2.601141e-09 |
| 16 | 0.08289699 | 0.01202087 | 0.0002411654 | -2.352252e-11 |
| 32 | 0.02852458 | 0.006354801 | 0.0001278255 | -1.385558e-13 |

This is evolution of an initialized equilibrium, not creation of an FRC. The finite-time drift and mesh/time controls do not establish long-time stability. Density is a mathematical choice, not a validated thermodynamic closure near the zero-pressure axis and surface.

## Time and prescribed-bias work

At nr=8, time-cap refinement uses identical spatial data. Material CFL can impose a smaller step; accepted step counts below show the actual amount of work. Every 0.004 s diagnostic boundary is reached exactly. The response residual is kinetic plus gas internal plus all-space induced magnetic energy, minus initial energy and integrated prescribed-source work; it is diagnosed without repair.

| Fixture | Time cap [s] | Steps | Source work [J] | Response residual [J] | Maximum speed [m/s] |
|---|---:|---:|---:|---:|---:|
| expansion | 0.002 | 40 | 0 | -1.9750221e-08 | 0.08360177 |
| expansion | 0.001 | 40 | 0 | -1.2154177e-08 | 0.0836018 |
| expansion | 0.0005 | 80 | 0 | -3.0385414e-09 | 0.08360184 |
| expansion | 0.00025 | 160 | 0 | -7.5963547e-10 | 0.08360185 |
| frc | 0.002 | 40 | 0 | -2.7466029e-09 | 0.03278401 |
| frc | 0.001 | 40 | 0 | -2.6011406e-09 | 0.03278417 |
| frc | 0.0005 | 80 | 0 | -6.5028161e-10 | 0.03278598 |
| frc | 0.00025 | 160 | 0 | -1.6257218e-10 | 0.03278643 |
| source | 0.002 | 30 | 1.9920148 | 0.0013077716 | 0.2908891 |
| source | 0.001 | 40 | 1.9928613 | 0.00044865824 | 0.2909376 |
| source | 0.0005 | 80 | 1.9931929 | 0.00011215771 | 0.2909577 |
| source | 0.00025 | 160 | 1.9932758 | 2.8039001e-05 | 0.2909627 |

The source fixture starts in the FRC equilibrium appropriate to B(0), then applies B(t)=sqrt(mu0) exp(-(t-0.04)²/(2*0.04²)) T through its peak at 0.04 s. The source is a prescribed spatially uniform magnetic bias, not arbitrary circular coils or a circuit model. The exact source-power formula is checked independently by differentiating magnetic energy at fixed material coordinates. Its driven flow has no analytic reference, so adjacent-grid differences remain necessary.

| Source radial cells | Max speed [m/s] | Max material displacement [m] | Source work [J] | Response residual [J] |
|---:|---:|---:|---:|---:|
| 4 | 0.1664504 | 0.002381477 | 2.0219054 | 0.0001139011 |
| 8 | 0.2909577 | 0.004109702 | 1.9931929 | 0.0001121577 |
| 16 | 0.5249064 | 0.007591669 | 1.9831922 | 9.792502e-05 |

### Saved-state comparisons

Differences use corresponding material nodes at matching physical saved times. The table reports the maximum across midpoint and final states, with the first row's mass weights; relative differences in tiny artificial equilibrium motion are not relative physical error.

| Test / reference | Velocity relative L2 | Velocity RMS difference [m/s] | Maximum position difference [m] |
|---|---:|---:|---:|
| source_4 / source_8 | 0.4377842 | 0.05215656 | 0.001728225 |
| source_8 / source_16 | 0.4716661 | 0.07293265 | 0.003481967 |
| frc_far_matching / frc_16 | 0.1220833 | 0.0003317907 | 3.407555e-05 |
| frc_modes4 / frc_16 | 0.000786327 | 2.143303e-06 | 2.359288e-07 |
| frc_quadrature7 / frc_16 | 6.007372e-10 | 1.647544e-12 | 3.758382e-13 |
| expansion_time_0.002 / expansion_time_0.001 | 3.256171e-07 | 2.017444e-08 | 6.672649e-10 |
| expansion_time_0.001 / expansion_time_0.0005 | 3.869079e-07 | 2.400807e-08 | 8.144916e-10 |
| expansion_time_0.0005 / expansion_time_0.00025 | 9.672667e-08 | 6.002005e-09 | 2.036222e-10 |
| frc_time_0.002 / frc_time_0.001 | 1.177595e-06 | 9.325748e-09 | 1.164514e-09 |
| frc_time_0.001 / frc_time_0.0005 | 1.31577e-05 | 1.08567e-07 | 2.76595e-08 |
| frc_time_0.0005 / frc_time_0.00025 | 3.288718e-06 | 2.713928e-08 | 6.914074e-09 |
| source_time_0.002 / source_time_0.001 | 0.0001706069 | 1.507816e-05 | 3.967816e-06 |
| source_time_0.001 / source_time_0.0005 | 7.028102e-05 | 6.21184e-06 | 1.553656e-06 |
| source_time_0.0005 / source_time_0.00025 | 1.756971e-05 | 1.552939e-06 | 3.883757e-07 |

The driven 8-to-16 mesh change is 47.167% in velocity L2. This is not a converged driven-flow result, even when source-work error is much smaller. Ideal surface-current response launches a narrow moving disturbance; spatial refinement and the later independently verified resistive extension must resolve it before interpreting formation. The current coarse-grid drive is a source/work and sensitivity fixture only.

## Flux inside material cells

The solver preserves psi at each moving non-axis vertex. It does not preserve every Q1-interpolated material-point psi exactly. For the independent prescribed deformation r_new=r(1+0.2z), z_new=z+0.1r² and initial a=0.5(1-0.2r²+0.1z) T, the actual interior discrepancy is:

| Radial cells | Material-point flux relative L2 drift |
|---:|---:|
| 4 | 0.0038648619 |
| 8 | 0.00092970117 |
| 16 | 0.00022992817 |
| 32 | 5.7322161e-05 |

This second-order spatial defect is numerical interpolation error, not physical resistivity or reconnection. For the evolved FRC/source cases:

| Case | Max vertex psi drift [Wb/radian] | Max material-point flux relative L2 drift |
|---|---:|---:|
| frc_4 | 5.421e-20 | 0.0006526302 |
| frc_8 | 5.421e-20 | 6.745071e-05 |
| frc_16 | 8.1315e-20 | 9.955069e-06 |
| frc_32 | 5.421e-20 | 1.852842e-06 |
| source_4 | 1.3553e-20 | 0.001310666 |
| source_8 | 2.7105e-20 | 0.0008936358 |
| source_16 | 2.7105e-20 | 0.0006319954 |

## Budgets, interpretation and remaining work

All sampled material mass changes are 0 relative. The maximum sampled vacuum relative residual is 2.897e-13. Minimum sampled density and pressure are 0.9628919 kg/m³ and 3.944536e-06 Pa, without floors. These sampled extrema do not assert bounds between outputs. Every accepted step validates positive geometry and finite state; rejected proposals retain the prior accepted state, time and work.

The evolving material interface is a new capability. The measured ideal interpolation defect, discrete equilibrium drift, spurious translation force and matching dependence remain explicit. Small budget and elliptic residuals do not prove complete physical accuracy. Q1/P0 material meshes can distort and have low-energy modes; there is no shock capture or remapping in this reference.

Next extend the compatible ideal force/vacuum treatment to physical resistive induction and positive Ohmic heating with independent diffusion/energy tests, then spatial coil sources and resolved driven-flow comparisons. Continue transport/ionization/Hall applicability, actual geometry/end conditions, shock/mesh robustness and converged formation diagnostics. An initialized reversed-field torus, a closed contour, numerical survival and particle confinement are distinct claims.

Reproduce measurements with `python -m examples.moving_static_measurements` and `python -m examples.moving_measurements --reuse`; assemble only complete source-matched records with `python -m examples.verify_moving`. The finest evolution cases can take several minutes. Sources for continuum concepts and the complete derivation are linked in MOVING_PLASMA.md; no published stability guarantee is applied to this discretization.

![Independent convergence](../outputs/moving_verification/convergence.png)

![Computed material fields](../outputs/moving_verification/material_fields.png)
