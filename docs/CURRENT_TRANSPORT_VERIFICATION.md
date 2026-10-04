# Current-transport verification — 2026-10-04

Frozen calculation source:
`2f682eece854d54b3b0d33026709f8d34617cc4633d8b3875a0f3e2b8e0d9540`.
All 16 cases completed and replayed arrays/metrics bitwise. Independent
root-polynomial/FFT integration at orders 16 and 20 measures actual continuous
fields, including force between points and the physical energy derivative.
The full suite passes **408 tests in 228.47 s**. All 13 older calculation source
closures remain unchanged. Evidence is in `outputs/current_transport_verification/`
and `outputs/current_transport_diagnostics/`; see the six-panel PNG/PDF report.

## Original local terms

The protocol varies fluid degree 3/5, radial elements 4/8/16/32, quadrature
8/12/16, field resolution separately, axial point count, density contrast, and
physical density 1e16/1e18/1e20 m^-3. Field polynomial degree is five. Matched
field resolutions have twice as many radial/axial field cells as fluid elements.
This is an instantaneous prescribed fixture in a filled periodic cylinder,
not FRC formation or complete species transport.

Degree-five relative L2 electric errors against the independent original terms:

| Fluid elements | Advection | Stretching | Self drift |
|---|---:|---:|---:|
| 4 | 7.17983e-2 | 6.35015e-2 | 2.99123e-2 |
| 8 | 5.26488e-3 | 9.54924e-3 | 1.88642e-3 |
| 16 | 2.34840e-4 | 4.73848e-4 | 7.73915e-5 |
| 32 | 8.45095e-6 | 1.72054e-5 | 2.71144e-6 |

Independent continuous drift-force component errors (not just nodal errors):

| Fluid elements | Radial | Azimuthal | Axial |
|---|---:|---:|---:|
| 4 | 4.30945e-4 | 6.04523e-3 | 3.53885e-2 |
| 8 | 1.27612e-5 | 4.25928e-4 | 2.47460e-3 |
| 16 | 5.00998e-7 | 2.74359e-5 | 1.59049e-4 |
| 32 | 2.55272e-8 | 1.72774e-6 | 1.00102e-5 |

At 32 elements, first-element component errors are 2.35505e-8, 1.62848e-6,
9.73869e-6; last-element errors are 1.07629e-6, 2.15170e-6, 6.36837e-6.
Both boundaries refine. Degree-three force converges more slowly: finest
combined global error 4.33069e-5 and last-element relative error 7.28561e-3.
Local relative norms near a vanishing wall force are more demanding; maximum
errors are recorded separately rather than concealed by a global norm.

Omitting interface terms retains the algebraic skew energy identity but makes
local forces inaccurate. At 32 degree-five elements the omitted-face advection
error is 1.32919e-3 versus 8.45095e-6 with faces. Separate field-grid controls
show substantial field-reconstruction error; changing fluid resolution alone
cannot eliminate it. No diffusion or energy repair was introduced.

## Continuity and energy measurements

At the finest degree-five resolution, strong continuity relative defect is
5.59788e-9, scaled div(J) is 4.40322e-7, and scaled axial normal-current jump is
1.37814e-7. Radial jumps are 7.36197e-15 and wall normal current 2.22446e-16;
wall radial velocity is zero. The nonzero defects are retained in interpretation
of the weak skew transport and must remain controlled in future coupled work.

Across all cases, advection/density-work pairing error is <=2.46009e-15;
stress/stretching pairing <=1.26520e-15; self-drift work backward error
<=3.85273e-16. Scaled axial force/actual angular torque are <=3.90582e-16 and
2.12385e-14. Original induction residual is <=1.61608e-11. Grid instantaneous
total energy backward error is <=3.04005e-14.

The independent physical current/electric/magnetic integrals agree with
production within 2.75850e-16 of the absolute total physical-power scale.
For the nonzero power comparisons, maximum relative discrepancy is
1.91576e-13. Order-16/20 power changes are <=9.57811e-17 on the same scale.
Relative errors of nearly zero self-drift work are not informative; absolute
and scaled differences are retained in the audit.

**Grid conservation is not exact continuous energy conservation.** The
independent degree-three physical total-power errors at 4/8/16/32 elements are
2.35289e-5, 3.81597e-7, 5.95719e-9, 3.72857e-11 W. Their largest backward error
is 7.34124e-9 and decreases with refinement. Degree-five physical residuals
are about 5.4e-11–9.6e-11 W, consistent with the original induction solve's
remaining residual. At the finest case the previously missing density exchange
-0.00149846144 W is paired by current advection; the full instantaneous physical
rate is -5.37703e-11 W. No heat was inserted to obtain that pairing.

## Scope and next gate

This completes frc-7mm.3 as an audited component. It does not replace an old
production solver, validate coupled time evolution, or finish the FRC goal.
Species pressure/current entropy and their common charge/particle-flux contract
are next (frc-7mm.4). Unforced complete-parent trajectories, heated moving exact
vacuum, physical closures/ionization, boundaries and formation remain open.
The pre-freeze audit sampler's periodic wrap test initially failed and was
corrected without changing production transport; the failure log is preserved.
