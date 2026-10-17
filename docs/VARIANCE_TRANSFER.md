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

## Mechanical resolution and actual transport heating

`outputs/reaction_thermal_transport` completes all three full actual 5 ns
transport stages (69881 EXIT0), starting from the captured original third-interval
first source. Only Ar8 velocity resolution changes; all old knots, density cuts,
populations and physical masses remain. These are stage calculations, deliberately
without the second reaction half, not new complete coupled 20 ns trajectories.

| Ar8 velocity cells | Final maximum T (eV) | Maximum excess over adiabatic T (eV) | Volume RMS excess (eV) |
|---|---:|---:|---:|
| 32 | 7.041390 | 6.726178 | 0.701676 |
| 44 | 0.390910 | 0.110259 | 0.012617 |
| 46 | 0.347248 | 0.061844 | 0.007697 |

Each trajectory completes 223 steps; relative transport-energy error is about
5.381e-11. The diagnostic compares with U*h^(gamma-1)=constant for reversible
ideal-gas compression. Numerical viscosity can represent shock entropy, so
excess heat is not by itself proof of an unphysical solution. The strong mesh
sensitivity is unresolved error evidence; these grids do not establish an
asymptotic order. Saved states, actual dt, per-cell fields and figures remain.

Explicit `velocity_refinement_rounds` enriches the mechanical mesh at existing
central-half density faces while retaining old knots. Capacity applies to whole
proposed rounds. Default zero is unchanged. Optional species selection allows
heavy-only enrichment; all other species still evolve. No density coarsening,
heat floor, population cutoff or changed mass is introduced.

Full 20 ns / 10 ns-source enriched run `reactive_material_argon_enriched_10ns`
is EXIT0: 1977 steps, 198.80 s, 21347 density cells, relative energy error
7.250e-11. Heavy meshes end at 46 cells (Ar9 24), electron at 50. Enriching the
electron mesh late reduces its acoustic bound and increases cost substantially.
The heavy-only variant `reactive_material_argon_heavy_enriched_10ns` is also
EXIT0 (confirmed by completed command record): 971 steps, 75.50 s, 21349 density
cells, energy error 2.842e-11. Its electron mesh remains 32. Relative to the
unenriched candidate its E/J differences are 5.271e-7 / 2.095e-5. Ar8 temperature
RMS/max differences are 0.008756 / 0.069926 eV; no local-convergence claim.
The heavy-only source SHA is
`9ea988bbfd2851b3ed8b828d5df837fffb4d1a5d0dade15c6f4b3befadc4d534`.

## Current discrepancy against the existing fine trajectory

`outputs/reaction_variance_fine_reference` reuses the already completed direct
512/1024-cell trajectories, retaining all fine cuts in union norms. Both use
physical argon/electron masses and the same 20 ns initial plasma profiles.
The reference omits reactions and the field-uncoupled neutral component; it is
not an exact or experimental solution, and its own uncertainty remains visible.

| Full endpoint | Relative E difference from 1024 | Relative J difference from 1024 |
|---|---:|---:|
| 32, variance | 10.390% | 66.305% |
| 64, variance | 4.993% | 57.312% |
| 128, constrained | 2.447% | 34.203% |
| 32 initial, all enriched | 10.615% | 65.774% |
| 32 initial, heavy enriched | 10.390% | 66.305% |
| 512, direct | 0.3303% | 8.0781% |

Late electron enrichment barely reduces its current discrepancy. The next
actual coupled experiment uses 128 electron cells from the start and 32 heavy
cells, with the same plasma-phase cap as the previous uniform128 case. The
new driver option records both resolutions; it does not change physical masses.
No claim that current is resolved or that the reference is converged.

## Preserved failures and arithmetic diagnosis

`precise_material_enriched` terminates EXIT1: enriched neutral decay rejects the
original non-subcycled step before any coupled time because source return
reduces the acoustic bound. Its charged affine control completes. The explicit
subcycled replacement retains the complete original durations and every actual
transport step; it is not a hidden retry or shortened successful run.

`reactive_material_argon_heavy_enriched_5ns` (22186 EXIT1) completes only 10 ns
of its requested 20 ns. The next first-source return creates six positive Ar2
amounts near 1.82e-318 in normalized coordinates. Their stored source heat is
exactly zero. Computing gamma*p*h/(m*N) gives 0/0 after m*N underflows, and the
acoustic eigenproblem rejects. `outputs/reaction_acoustic_underflow` reconstructs
the actual source and verifies the zero source heat independently with 100-digit
decimal integration. Reordering sound to gamma*(gamma-1)*(U/N)/m returns the
correct zero without adding heat. Temperature and numerical-viscosity ratios
need the same careful ordering. Finite-precision source chemistry at such tiny
populations is still a limitation; zero stored heat is not a physical accuracy
claim. The original failed source, knots and accepted prefix remain unchanged.

## Live calculations and next action

At the current continuation, poll the actual handles below; older live entries
elsewhere are historical. No result is complete until the full target is reached.

- **73344**: `reactive_material_argon_variance_5ns`, frozen
  `/private/tmp/frc-variance-20261017`, 15 ns accepted; fourth interval active.
- **57807**: `precise_material_variance`, same frozen root; neutral16/32 complete,
  neutral64 active, followed by three charged affine cases.
- **47880**: `precise_material_enriched_subcycled`, frozen
  `/private/tmp/frc-enriched-subcycled-20261017`; full neutral decay active,
  then full charged affine control. Initial16, one enrichment round, capacity64.
- **77325**: `reactive_material_argon_e128_h32_10ns`, main checkout source;
  10 ns accepted, second interval active. Do not edit its imported source until
  terminal. Compare its full endpoint with the existing uniform128 and fine
  reference to evaluate independent species resolution.
- **16195**: `reactive_material_argon_factored_heavy_5ns`, frozen
  `/private/tmp/frc-acoustic-factored-20261017`; complete20ns request with
  heavy-only enrichment and algebraically reordered sound/temperature/viscosity.
  This prototype is not yet copied into main `precise_material.py`.

Keep the old frozen roots intact. Compare complete 10/5 ns local temperature,
velocity, fronts and E/J when available. No 2.5 ns run yet: resolve the present
thermal/arithmetic errors before further mesh growth. The 20 ns planar component
is not the 1 m lab. Gas chemistry/collisions, heat, driven coupling, boundaries
and independent physical formation evidence all remain required.
