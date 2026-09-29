# Electron-current inertia and physical regime checks

The new `MaterialCurrent` component implements fixed-material toroidal induction with physical electron inertia, scalar prescribed relaxation and genuine exterior magnetic vacuum. It separately accounts for induced magnetic energy, electron drift kinetic energy, positive collisional heat and prescribed-source work. It is a linear component reference, not a moving Hall/thermal solver or an FRC formation result. The prior heated-edge incompatibility remains open. See [the equations and discrete energy derivation](CURRENT_INERTIA.md) and [conditional transport conventions](TRANSPORT_REGIME.md).

All 35 declared component cases are complete. Source/configuration and lossless array hashes are checked; report, time-audit, scan and rendering sources are separately recorded. Radius is 1 m, analytic field coefficient/bias is 1e-8 T, and the constant-profile decay uses D=0.1 m²/s. Density follows physical de²=me/(mu0 ne e²). These are explicit component inputs, not selected laboratory operating parameters.

The de=0 rows explicitly select the massless comparison; no finite carrier density is assigned to them.

## Independent continuum comparisons

Sphere decay has continuous dipole/vacuum matching and exact rate D k²/(1+de² k²). Collisionless screening has an independent modified-spherical-Bessel solution and vanishes continuously with carrier density. A third exact sphere correlates prescribed density/resistivity profiles so the same continuum mode has a known decay; this verifies variable coefficients without a numerical forcing or an energy-only assertion. It is a manufactured profile, not a species closure.

| Case | Magnetic field relative L2 | Current relative L2 | Local heat relative L2 | Initial current relative L2 |
|---|---:|---:|---:|---:|
| decay_d0_8 | 0.04302239 | 0.01317183 | 0.0297218 | 0.02880242 |
| decay_d0_16 | 0.02107254 | 0.00330263 | 0.007356365 | 0.008975305 |
| decay_d0_32 | 0.0104809 | 0.0008264379 | 0.001835126 | 0.002923812 |
| decay_d0.05_8 | 0.04293943 | 0.01383523 | 0.03076228 | 0.02880242 |
| decay_d0.05_16 | 0.02106313 | 0.003468075 | 0.007797864 | 0.008975305 |
| decay_d0.05_32 | 0.01047975 | 0.0008686021 | 0.002025159 | 0.002923812 |
| decay_d0.2_8 | 0.04253866 | 0.02497522 | 0.03645654 | 0.02880242 |
| decay_d0.2_16 | 0.02101163 | 0.007687326 | 0.01044695 | 0.008975305 |
| decay_d0.2_32 | 0.01047327 | 0.002479929 | 0.003220468 | 0.002923812 |
| decay_d1_8 | 0.04243365 | 0.02877453 | 0.0385463 | 0.02880242 |
| decay_d1_16 | 0.0209984 | 0.008966779 | 0.01112837 | 0.008975305 |
| decay_d1_32 | 0.01047161 | 0.002921072 | 0.00345751 | 0.002923812 |
| decay_d0.2_64 | 0.005232557 | 0.0008288658 | 0.001048592 | 0.00098355 |
| screen_d0.05_8 | 0.4623484 | 0.2285301 | n/a (eta=0) | zero initial state |
| screen_d0.05_16 | 0.2735415 | 0.07059533 | n/a (eta=0) | zero initial state |
| screen_d0.05_32 | 0.1441267 | 0.01888342 | n/a (eta=0) | zero initial state |
| screen_d0.2_8 | 0.07945078 | 0.02772173 | n/a (eta=0) | zero initial state |
| screen_d0.2_16 | 0.0399022 | 0.006986867 | n/a (eta=0) | zero initial state |
| screen_d0.2_32 | 0.01998004 | 0.001750552 | n/a (eta=0) | zero initial state |
| screen_d1_8 | 0.004171138 | 0.002054792 | n/a (eta=0) | zero initial state |
| screen_d1_16 | 0.001934158 | 0.0005188341 | n/a (eta=0) | zero initial state |
| screen_d1_32 | 0.0009467435 | 0.0001300547 | n/a (eta=0) | zero initial state |
| screen_d4_8 | 0.0002667683 | 0.000136995 | n/a (eta=0) | zero initial state |
| screen_d4_16 | 0.0001224482 | 3.462788e-05 | n/a (eta=0) | zero initial state |
| screen_d4_32 | 5.974016e-05 | 8.682444e-06 | n/a (eta=0) | zero initial state |
| screen_d0.05_64 | 0.07307775 | 0.004807223 | n/a (eta=0) | zero initial state |
| profiled_8 | 0.04245396 | 0.0279506 | 0.04451081 | 0.02880242 |
| profiled_16 | 0.02100089 | 0.008719013 | 0.01310823 | 0.008975305 |
| profiled_32 | 0.01047192 | 0.002843038 | 0.004149714 | 0.002923812 |
| decay_d0.2_16_time | 0.02101163 | 0.007687331 | 0.01044688 | 0.008975305 |
| decay_d0.2_16_matching | 0.02101278 | 0.00807018 | 0.01099183 | 0.009429044 |
| decay_d0.2_16_angular | 0.01190017 | 0.006597619 | 0.009068567 | 0.007697726 |
| decay_d0.2_16_quadrature | 0.02101163 | 0.007687326 | 0.01044695 | 0.008975305 |
| decay_d0.2_16_grading | 0.02108228 | 0.02828312 | 0.04767361 | 0.033769 |
| decay_d0.2_16_modes | 0.02101163 | 0.007687375 | 0.01044695 | 0.008975365 |

The Q1-gradient magnetic field converges approximately at first order. Inertia retains weak-current error from the interpolated initial a, especially near the surface; at de/R=0.2 the nr=8/16/32/64 current errors are about 0.02498/0.007687/0.002480/0.0008289. This approaches roughly order 1.5, rather than the near-second-order current error after strong high-mode damping in the massless comparison. The initialized field/current inconsistency is measured at t=0 and not removed. Surface grading alone worsens this current error in the recorded control.

## Thin screening and false reversal

The exact collisionless sphere has positive interior Bz everywhere for positive applied bias. Its numerical minimum and negative-field volume are therefore direct artifact checks; a small energy residual must not turn them into an FRC detector.

| Case | Minimum Bz / bias | Negative-Bz material volume [%] | Largest unphysical positive Jtheta [A/m²] | Central Bz / bias |
|---|---:|---:|---:|---:|
| screen_d0.05_8 | -0.003303198 | 40.26677 | 0.000978341 | -4.950709e-11 |
| screen_d0.05_16 | 1.717459e-08 | 0 | 0 | 1.692482e-08 |
| screen_d0.05_32 | 5.744571e-08 | 0 | 0 | 5.726324e-08 |
| screen_d0.2_8 | 0.06193851 | 0 | 0 | 0.06175256 |
| screen_d0.2_16 | 0.06571295 | 0 | 0 | 0.06566316 |
| screen_d0.2_32 | 0.06688808 | 0 | 0 | 0.06687535 |
| screen_d1_8 | 0.850256 | 0 | 0 | 0.8501584 |
| screen_d1_16 | 0.8506016 | 0 | 0 | 0.8505761 |
| screen_d1_32 | 0.8508007 | 0 | 0 | 0.8507943 |
| screen_d4_8 | 0.989624 | 0 | 0 | 0.9896169 |
| screen_d4_16 | 0.9896391 | 0 | 0 | 0.9896373 |
| screen_d4_32 | 0.9896511 | 0 | 0 | 0.9896506 |
| screen_d0.05_64 | 7.518742e-08 | 0 | 0 | 7.512947e-08 |

The unresolved de/R=0.05, nr=8 mesh produces false reversal and opposite-sign current. Both disappear in the recorded finer samples, but thin-layer magnetic-field error remains appreciable. Even positive field samples do not by themselves prove topology accuracy between samples. No clipping, density floor, resistance enhancement or heat cutoff hides the coarse artifact.

## Time, coefficients and boundary controls

| Same material mesh: baseline / reference | Current difference [%] | Cell heat-density difference [%] |
|---|---:|---:|
| decay_d0.2_16 / decay_d0.2_16_time | 9.262407e-07 | 1.173497e-05 |
| decay_d0.2_16 / decay_d0.2_16_matching | 0.04539906 | 0.06053161 |
| decay_d0.2_16 / decay_d0.2_16_quadrature | 4.221485e-10 | 5.713064e-10 |
| decay_d0.2_16 / decay_d0.2_16_modes | 3.25707e-05 | 2.780887e-05 |

Angular/grading cases have different material meshes and remain separate continuum-error rows above; no equal-index field comparison is made across them. The explicit fixed-spatial-eigenmode audit isolates temporal accuracy independently of continuum spatial error:

| dt [s] | Exact semidiscrete current relative L2 error |
|---:|---:|
| 0.04 | 0.00016996475 |
| 0.02 | 4.2467605e-05 |
| 0.01 | 1.0615429e-05 |

Time-error halving ratios are 4.002221 and 4.000555. This confirms second order on this resolved fixed mode; midpoint remains unsuitable for arbitrarily underresolved dissipative steps.

## Energy and current response

Independent true-sphere volume integrals of B²/(2 mu0), me J²/(2 ne e²) and local heat include the exact infinite exterior dipole energy. Orders 48 and 96 agree, and uniform-decay totals independently match their closed formulas. Nonuniform-coefficient and deformed-geometry tests also evaluate physical energies separately from the mixed matrices.

Maximum continuum energy-quadrature relative change is 9.172e-14. Maximum whole-run response residual divided by the corresponding continuum magnetic/current/heat energy scale is 8.414e-15; maximum equation residual is 2.341e-12. Positive heat is retained per cell. These budget checks do not erase the measured field/current errors or establish a physical surface closure.

## Conditional argon and local-edge assessment

The original zero-current Ar+ reference (ne=1e20 m^-3, Te=5 eV, Ti=1 eV, radius=0.035 m, B=0.005 T) has di/R=4.0972, de/R=0.015183, and conditional electron mean-free-path/R=0.11088. Thus the supplied global parameters do not justify Hall-free MHD. Smaller edge-gradient lengths are more restrictive.

Its conditional Coulomb momentum resistivity is 8.575415e-05 ohm m; the separate Z=1 parallel/unmagnetized estimate is 4.375212e-05 ohm m. Prescribed conductivity is unchanged. The two coefficients have distinct limits; neither replaces the full tensor, neutral collisions or thermal transport.

The archived scan contains 43 explicitly assumed states and 41 samples of a prescribed vanishing-density edge. Its fixed current is a diagnostic scale, not a self-consistent solved edge. It exposes rising Hall/inertial/Debye and drift ratios as vacuum is approached. Fixed charge/ionization, Maxwellian distributions and temperatures are assumptions, not predictions; no tested state is labeled fully physically validated.

Next integrate the required Hall/toroidal and electron/ion thermal equations with a justified edge treatment, preserving genuine exterior vacuum and independent wave/current/energy checks. The linear inertia component alone must not be inserted into moving dynamics without its geometric current-energy force and convective terms. Spatial coil work, actual gas ends, robust meshes/shocks and converged formation remain required. Stage 2 remains gated.

Reproduce with `python -m examples.current_measurements --reuse`, `python -m examples.transport_study`, then `python -m examples.verify_current`. The numerical suite currently has 192 passing tests. All evidence is in `outputs/current_verification`.

![Current accuracy](../outputs/current_verification/current_accuracy.png)

![Conditional transport scales](../outputs/current_verification/transport_scales.png)
