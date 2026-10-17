# Conservative thermal-transfer candidate, 2026-10-17

The lab-facing simulator goal is active and the physical FRC goal remains open.
This continuation targets the local thermal error in bead 5.5.5. The new
`variance` option is experimental; Galerkin/local_total remains the default.
Read `plasma_bench/variance_projection.md` for the equations and contracts.

The candidate fits velocity in occupied volume with exact block momentum. It
selects between two energy-conserving scalings of that shape by volume error,
then integrates source heat plus the nonnegative squared velocity discrepancy.
This avoids recovering target heat by subtracting two mismatched kinetic
energies. Local total energy is redistributed with momentum and measured.
Numerical mixing is not physical reaction heating. A constant selected shape
can erase velocity structure; positivity/conservation alone are insufficient.

## Actual source and complete-run evidence

- `outputs/reaction_thermal_source`: complete 10–15 ns replay from the original
  accepted 10 ns state, 223 actual steps, 312.918 s, EXIT0. Every saved endpoint
  array equals the original 15 ns state. Both evolved quadrature sources and
  old velocity knots are saved. Replay source SHA
  `ec5eeb9730cd7234589ded15e4d3209cb626f304ed054c202f426af4eb939504`.
- Its second source already has Ar8+ temperature 0.165279–7.041390 eV. The old
  constrained return produces 2.22e-12–14.221326 eV; source-to-target T RMS is
  0.366200 eV, maximum change 7.182864 eV. The variance candidate gives
  0.186225–8.381399 eV, RMS 0.054561 eV, maximum change 1.342937 eV. Velocity
  RMS decreases 198.98 to 181.08 m/s. This isolates one transfer and cannot
  remove heat already present in its input. All narrow cells/species remain.
- `reactive_material_argon_variance_10ns`: full 20 ns, EXIT0, 971 actual steps,
  21351 final density cells, 74.349 s. Relative energy error 2.84206e-11;
  largest final-return T change 0.046388 eV. Final trace-ion minima remain
  around 0.183 eV (Ar9+ 0.19856 eV). Source SHA
  `501fb7136404f9a442241b75239c725cbda592ac0273607bfa52166ba17ff322`.
- `reaction_variance_comparison`: against the completed matching original
  10 ns-source run, full-endpoint E/J relative differences are
  1.19020e-6 / 3.50460e-4. Per-species local norms/front shifts are retained.
- `reaction_remap_variance`: all full 20 ns zero-reaction controls complete,
  EXIT0, 56.00 s. Direct/10 ns/5 ns steps are 970/971/972. E differences from
  direct transport are 3.14214e-7 / 3.09182e-7; J differences are
  1.02096e-4 / 2.57461e-4. Observed mixing heat is zero at saved precision;
  redistribution L1 is 2.806e-22 / 6.183e-22 J/m².
- `reactive_material_argon_variance_n64_10ns`: full 20 ns, 71873 EXIT0,
  2041 actual steps, 42547 density cells, 276.423 s, relative energy error
  5.55259e-12. Maximum final-return T change is 0.071675 eV. Against the
  32-cell variance run, E/J differences remain 6.65162% / 32.21364%; Ar8+
  T RMS/max differences are 0.034847 / 0.502636 eV. Electron-current spatial
  accuracy is still unresolved. The plasma-phase cap also halves, as in the
  original driver, so this is not isolated spatial-order evidence.

## Independent limits and unresolved local accuracy

`reaction_variance_analytic` uses independently integrated 100-digit moments
for the density-contrast counterexample. The nonzero root would recreate
velocity endpoints -0.25 and 14999.75 m/s; comparing both roots selects the
correctly derived constant 0.5 m/s alternative. The code matches it exactly
at float precision. This selection deliberately loses unresolved structure.

For uniform density and u=1000*x², 8/16/32/64 elements give velocity L2 errors
4.65847/1.16462/0.291155/0.072789 m/s: observed order 2. Temperature mixing
errors have observed order 4; finest maximum change is 7.3121e-10 eV. Independent
adaptive integration of the continuous squared velocity error agrees with the
mixing heat. Relative total-energy errors are below 3e-17. These are mathematical
transfers, not physical argon runs or a general positivity/accuracy theorem.

`reaction_variance_refinement` retains the earlier rejected source and all old
knots. Neutral maximum T changes decrease overall 0.02697→0.00527 eV and
volume RMS 0.001244→0.000250 eV. But its velocity error is larger than the old
constrained projection; trace Ar9+ maximum thermal error is not monotone.
Nearest-face targets in this MPFR script are evaluated outside the high-precision
wrapper, affecting which nearly tied faces are added; actual counts are recorded
and no original density face is deleted. No asymptotic real-source claim.

## Live calculations and next action

At 02:55 UTC, full 5 ns-source 20 ns run **73344** is live in
`outputs/reactive_material_argon_variance_5ns`, accepted through 15 ns and
starting the fourth interval. The third interval completes in 309.21 s with
223 steps, 227825 density cells and relative energy error 2.677e-11. Its Ar8+
temperature range is 0.18625–8.66800 eV, versus the original
2.22e-12–14.22133 eV. Largest last-return temperature change is 1.40180 eV.
This is improvement, not local thermal convergence. Source SHA
`eb54713073c61523dd878b035ea91f406ba57ced373f0409b95e1c3ae060c28f`.
The six full neutral/charged analytic expansions **57807** are live in
`outputs/precise_material_variance`; neutral16 is complete, neutral32 running.
The new 64-cell full candidate **71873** and comparison **96533** are terminal
EXIT0, with results above. This checks spatial sensitivity of the changed method.
All execute frozen `/private/tmp/frc-variance-20261017`; do not edit that root.
The 5 ns source differs from the earlier complete 10 ns candidate only by the
explicit exactly-zero-energy residual case, which does not affect those warm
trajectories. Each run archives its actual source independently.

Focused transport diagnosis **69881** is also live, driver
`outputs/reaction_thermal_transport/stage.py`, using the same frozen source.
It evolves the actual full 5 ns transport stage after the captured original
first source, then repeats with only Ar8 mechanical elements refined once/twice.
Every old knot and density cut stays. Compare final heat with the independent
ideal-gas invariant U*h^(gamma-1) to distinguish reversible compression from
numerical entropy production. These are explicitly transport-stage calculations,
not new coupled 20 ns trajectories. Preserve all outputs and actual exits.

Poll those handles to terminal completion, then compare full 10/5 ns local
temperature/velocity/front/E/J results. Preserve failures and actual clock.
No 2.5 ns run yet: fix remaining local accuracy before further mesh growth.
All original transport, gas chemistry, driven coupling, wall and experimental
validation requirements remain. This 20 ns planar component is not the 1 m lab.
