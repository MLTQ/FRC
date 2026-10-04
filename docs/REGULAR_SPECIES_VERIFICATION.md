# Regular species coupling verification — 2026-10-04

The final suite passes **416 tests in 210.84 s**. Eight focused species tests
pass in 7.84 s. All **24 declared cases** completed, replayed every saved array
and numerical metric bitwise on the original platform, and passed independent
physical reconstruction at quadrature orders 16 and 20. No study remains live.

Frozen calculation source:
`45a71e5b59cbcab6bb8556d531510d68da8389b16abcc82e1e61c0bc3237b78d`
(64 Python files, code/companion archive and runtime versions). Evidence lives
in `outputs/regular_species_verification/`; logs, source checks, PNG/PDF figure
and machine-readable summary are in `outputs/regular_species_diagnostics/`.
The independent local enthalpy audit has source
`4c6f66538180f6ef8a72d72726a867d5c06751b64eb25daeb60494f794d9281c`
and its own archive and results in `outputs/regular_species_enthalpy/`.

## Protocol and original-equation accuracy

Four meshes at fluid degree 3 and 5, quadrature orders 8/12/16, densities
1e16/1e18/1e20 m^-3, a density-contrast control and a finer field-grid control
give 14 instantaneous cases. Both temperatures vary independently. The field
degree remains 5; field cells subdivide fluid elements. The analytic oracle
differentiates the original pressure directly and uses the original species
advection equation, not the discrete work identity.

Independent continuous relative L2 errors at degree 5:

| Fluid elements | Electron entropy rate | Ion entropy rate | Electron pressure E | Ion pressure E |
|---:|---:|---:|---:|---:|
| 4 | 2.41562e-2 | 2.27232e-2 | 1.43126e-2 | 1.73734e-2 |
| 8 | 1.15355e-3 | 1.09068e-3 | 7.69591e-4 | 8.54753e-4 |
| 16 | 4.04405e-5 | 3.82809e-5 | 2.95661e-5 | 3.16903e-5 |
| 32 | 1.30090e-6 | 1.23178e-6 | 9.79306e-7 | 1.04547e-6 |

Finest electron first/last-element entropy errors are 1.31113e-6/1.12230e-6;
maximum error normalized by target RMS is 8.74405e-6. Corresponding ion values
are 1.31102e-6/1.02091e-6/8.24394e-6. Both species' pressure axis/wall/max
errors also decrease. Degree-3 finest entropy errors are 9.95929e-6 and
6.56150e-6; the lower-degree fluid representation remains distinguishable.
Changing audit quadrature changes integrated error metrics by at most
1.13970e-8 relative to those errors; sampled maxima use different locations
and are not included in that quadrature comparison.

## Conservation, physical energy and local enthalpy

Maximum species grid pressure/thermal work backward error is 1.32826e-16.
The pressure-induced magnetic plus current-inertia response gives maximum
grid total-power backward error 2.45753e-12; maximum checked induction
residual is 1.24843e-11 at tolerance 1e-10. The reconstructed full original
Ohm equation has maximum backward error 2.34994e-16. Uniform but unequal
temperatures at varying density exercise the analytic pressure curl kernel.

Grid conservation does not imply exact continuous thermodynamic conservation.
The degree-3 electron continuous thermal-minus-pressure-work discrepancy is
-6.95401e-9, -1.13341e-10, -1.78894e-12 and -2.73774e-14 W under refinement.
At degree 5 the finest discrepancy is -2.15106e-16 W in this smooth fixture.
Integrated entropy rates are normalized against their actual rate scales;
large absolute cancellation residues are not presented as relative errors.

An independently differentiated physical audit measures
`epsilon_t - J dot E_pressure + div(a*gamma*h*J/e)` within each polynomial
piece. Electron local relative L2 defect decreases
2.04643e-2 → 1.11992e-3 → 4.75271e-5 → 2.12690e-6; ion defect decreases
2.45076e-2 → 1.30586e-3 → 5.89135e-5 → 2.98064e-6. Finest electron/ion
maximum defect over term RMS is 1.38802e-5/1.88519e-5. The actual reconstructed
strong div(J) is included and its scaled norm is still 4.40322e-7 at the
finest mesh. Graph charge continuity does not make it identically zero.

The finest integrated piecewise enthalpy divergence is -3.73488e-7 W for
electrons and 1.06834e-11 W for ions. This is distinct from the tiny global
pressure/thermal work discrepancy. Distributional interface contributions
are not automatically included by integrating a piecewise strong derivative.
The full physical interface enthalpy/Poynting balance remains an adoption gate.

## Frozen-current evolution and stability

Eight unforced entropy-advection trajectories keep density and current
**externally fixed**. The exact solution is sheared periodic translation under
Jz=2G/mu0*(1-2r²/R²), with the physical species velocity a_s*Jz/(e*n).
The electron fastest-axis transit is 1.208029311 microseconds. The physical ion
fraction gives 0.08796960 s for the corresponding ion test; the ion mass is
not modified to shorten its transit. This is a verification interval, not a
prediction of a collisionless laboratory lifetime.

| Fluid elements | Continuous entropy-signal error | First-element error | Continuous quadratic relative drift |
|---:|---:|---:|---:|
| 4 | 5.78249e-2 | 7.11412e-3 | -2.95681e-7 |
| 8 | 1.87973e-3 | 5.31807e-6 | -3.82218e-10 |
| 16 | 3.05425e-5 | 4.26861e-8 | -3.51823e-13 |

Finest wall error is 6.17115e-5 and maximum error over signal RMS is
2.75903e-4. Initial continuous interpolation errors are below 6.43e-15 on
these meshes. Electron/ion and reversed-current controls retain the same
matched dimensionless transport error. The 120/240/480/960-step comparisons
give successive endpoint-difference ratios **15.99985 and 15.99998**, resolving
RK4 temporal error separately from the spatial error.

Finest grid quadratic drift is -1.78565e-15, distinct from continuous drift.
Finest physical entropy drift per particle is -5.06244e-17. Nonlinear internal
energy is not the preserved quadratic invariant: its relative drift changes
sign and decreases from -1.15760e-7 at 4 elements to 3.03419e-13 at 16 elements.
There is no thermodynamic energy repair.

Two complete small frozen spectra and randomized discontinuous-state controls
retain the failures that motivated the projection:

- The original point-transfer unskew operator has positive growth; at four
  elements max(Re lambda)/spectral radius is 0.0177381. Its exact exponential
  increases the selected state's quadratic norm by **22.1418%**.
- Naive point-transfer skew symmetrization removes that growth but retains
  order-unity local entropy-rate error under refinement.
- Physical-projection skew transport has normalized maximum real eigenvalue
  1.52282e-16, symmetric-part norm 4.00389e-16 and exponential quadratic drift
  -3.47709e-16. Dimensionless constant-state residual is 3.59994e-17.

Both species conserve entropy and pair thermal work in the discontinuous
control (work backward errors below 2.05e-16 in the four-element case). This is
not a discontinuous shock accuracy test or a nonlinear entropy inequality.

## Preservation and remaining work

All 15 earlier calculation source closures remain unchanged. The initial
24-case archive is retained in `outputs/regular_species_before_dtype_guard/`.
A final float64 input guard prevents integer/float32 rate truncation; all
24 arrays and numerical metrics remained bitwise identical after that change,
and the final 416-test suite includes its regression. The report was rendered
and visually inspected. No external job was launched, remote modified or
research data published.

This completes a verified **species component**. `frc-7mm.4` remains open for
combining actual bulk continuity/momentum, current advection/drift stress,
species pressure/entropy and magnetic evolution in unforced parent wave/contact
trajectories, including local interface/Poynting work, temporal and spatial
controls, and atomic failure. The original fluid order/spectral limitations,
heated exact-vacuum interface, selected physical closures and real coil/end
boundaries remain open before converged formation or experiment recommendations.
The user's 1 m quartz-tube envelope is recorded in RESEARCH_QUESTIONS.md;
the smaller verification fixtures do not select a laboratory configuration.
