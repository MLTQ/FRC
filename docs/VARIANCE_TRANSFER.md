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
| 128 electron / 32 heavy, variance | 2.4475% | 34.2028% |
| 512, direct | 0.3303% | 8.0781% |

Late electron enrichment barely reduces its current discrepancy. The complete
128-electron/32-heavy case77325 EXIT0 reaches20ns in4224steps,798.04s and84091
final density cells, relative energy error1.3853e-12. It uses the same common
plasma-phase cap as the earlier128-cells-per-species case. E/J differences between them
are0.036045%/0.056490%; the older case uses constrained return and the new one
variance, so heavy-species resolution is not the only changed choice. Observed
runtime is about25% lower than the prior1069.60s equal-count case, not a controlled
hardware benchmark. Its source SHA is
`10621ccac12ba7035ed9b69d3b3f0ccf5f55d5b1a427af961137dc1abf030ede`.

The small128-per-species versus e128/h32 discrepancy supports retaining the cheaper
heavy mesh for the256-electron run. Current still differs34.2% from the
existing1024-cell reference; that reference is itself unconverged. The updated
scientific comparison figure is inspected. Initial per-species resolution changes
neither physical masses nor support/vacuum rules.

Run31198 is now EXIT0 at full20ns in8647steps (4869/3778),3028.127s,
166537density cells, minimum width8.11966e-24m and Erel3.65072e-13. Initial
and final electron velocity counts are256. Fine comparison8942 EXIT0 gives
E/J differences1.02547%/19.07075% from1024, down from2.44748%/34.20280% at
128electrons. Electron-front difference falls3.39308→1.85810mm; the32-cell
heavy-ion front difference remains34.3293um. The existing512/1024 direct
comparison retains E/J0.33030%/8.07805%, electron-front0.771314mm, so the
reference uncertainty is not zero. No isolated spatial order: phase cap also
halves, and the new case includes separately measured arithmetic/scaling changes.
Scientific PNG and rendered PDF inspected. Thermal/current accuracy stays open.

An additional explicit experiment changes128 initial electron cells from cosine
to uniform spacing, retaining32 cosine heavy cells, variance return,128-bit
coordinate scaling and the same common phase cap. The continuous initial
density/temperature/support are unchanged; each cell contains n times its actual
width. Centre spacing changes245.412→156.25um, edge spacing3.012→156.25um.
This moves resolution toward the measured interior current discrepancy, but
can degrade the vacuum rarefaction. No evolved knot or fine density cut is
removed. Earlier "uniform128" meant equal cell counts with cosine spacing,
not this uniform-spaced grid. Full20ns completion and endpoint comparison are
required before any accuracy/cost conclusion.

The first attempt9098 now completes20ns in4027steps (2037/1990),1021.890s,
110129density cells and Erel−4.62663e-12. Its initial snapshot exposes a setup
error:128uniform density cells were given to the constructor without explicit
velocity knots; cosine-target selection chose only102velocity elements, with
widths156.25–312.5um. Thus this run never supplied the intended uniform128
mechanical mesh. Existing fine comparison13317 EXIT0 measures E/J differences
5.25086%/48.44737% and electron-front difference11.25444mm from1024, versus
2.44748%/34.20280%/3.39308mm for the cosine128-electron/32-heavy case.
This rejects the mixed102/128 mesh, not the intended uniform128 discretization.

The driver now supplies every initial electron face through `velocity_faces`
for the explicit uniform option and records actual initial velocity counts.
Fresh run91194 starts with[32,32,128] actual velocity cells, independently read
from saved indices; all128electron widths are156.25um within float spacing.
It repeats the full20ns request and preserves the original mixed-mesh output.
The comparison now reads actual initial density/velocity counts rather than
assuming the requested counts were realized. Cosine defaults remain unchanged.

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
decimal integration. Reordering sound to gamma*(gamma-1)*(U/N)/m evaluates the
stored cold state without 0/0. This alone is insufficient: the independent
full transfer integral below reveals mixing heat lost before sound evaluation.
Temperature and numerical-viscosity ratios need the same careful ordering. Finite-precision source chemistry at such tiny
populations is still a limitation; zero stored heat is not a physical accuracy
claim. The original failed source, knots and accepted prefix remain unchanged.

## Retained heat and binary coordinate scaling

The six cold source cells have nonzero projected velocity near405m/s. An
independent100digit integration of the declared variance-transfer formula gives
returned temperatures0.0226328089–0.0226687236eV, while the old stored U is zero.
Thus merely returning zero sound avoids the crash but retains a local thermal
error. `outputs/reaction_thermal_scaled` applies an explicit128-bit coordinate
shift: multiply stored relative N/momentum/U by2**128 and divide their common
physical population scale by2**128 before SI constants and tiny cell weights.
The physical amounts and equations do not change. No floor or population cutoff.

All six returned temperatures then agree with the independent decimal moments
to at most7e-18eV; normalized cell amounts become~6.18e-280 instead of~1.82e-318.
The actual admissible step is4.574096671e-11s, with no nonfinite sound values.
The option rejects an unrepresentable scale, retains default shift0, and cannot
make finite-precision source chemistry accurate at arbitrarily tiny populations.
Exact source identities are saved in provenance JSON. All archived inputs and
old failed results remain.

`reactive_material_argon_scaled_heavy_10ns` completes the full20ns (5892 EXIT0):
971steps,21349density cells, Erel2.842085e-11. Compared with the original
heavy-enriched10ns case, E/J differences are9.43e-15 /1.658e-8. This includes
both coordinate scaling and equivalent sound/temperature/Q arithmetic changes,
so it is not an isolated scaling study. Source SHA
`9bf93c336c75d48834f8ed91a9171f467919d951408f1fef11c6441b603ae467`.
The full5ns-source counterpart61204 is live from the same frozen source.
After77325 terminated, the exercised arithmetic and optional coordinate shift
were integrated into main with their companion contracts. Shift0 and default
projection choices remain unchanged; no broad formation adoption is claimed.
The unscaled factored-ratio case16195 now accepts15ns, with Ar8 maximumT0.65149eV
versus8.66800eV in the unrefined variance case. Neither prefix is a complete20ns
result or a local thermal-convergence claim.

## Completed independent expansions

`precise_material_variance` (57807 EXIT0,2619.44s) completes all six original
neutral-decay/charged-affine durations at16/32/64 initial velocity cells. Finest
neutral front error8.60e-16m, energy error5.78e-14 and population error3.44e-14;
finest charged front error6.94e-17m and energy error3.27e-15. The neutral exact
population law and independent affine trajectory are preserved. These fixtures
explicitly use8ULP cleanup; reacting argon uses zero. Figure inspected.

`precise_material_enriched_subcycled` (47880 EXIT0) completes the full6.083704us
neutral decay with987coupled intervals/1446actual transport steps, front error
1.16e-15m, energy error7.18e-14, population error4.33e-14. Its full4.446830us
charged control also completes75steps; no source/enrichment applies in that
control. Numerical thermal change is a signed-3.19e-17J/m2 roundoff-sized
observation, not negative physical reaction heating. The failed non-subcycled
attempt remains retained. These affine fixtures do not establish accuracy on
the nonlinear reacting trace-ion profiles.

## Live calculations and next action

The original variance5ns run73344 is now terminal EXIT0: full20ns,972 transport
steps (325/203/223/221),3818.456s,2415015 final density cells, minimum width
1.66804e-28m and Erel2.76935e-11. Source SHA
`eb54713073c61523dd878b035ea91f406ba57ced373f0409b95e1c3ae060c28f`.
The existing complete-endpoint comparison33358 EXIT0 now includes its entire
union. Against variance10ns, E/J differences are4.53748e-6/4.33645e-4. Neutral
front shift1.68884um, velocity RMS0.19978m/s, temperature RMS3.84577e-5eV;
electron front shift7.91207e-11m, velocity RMS0.33377m/s and temperature RMS
0.000208194eV. Trace Ar8 remains inaccurate: T RMS0.161226eV, maximum difference
2.231827eV, final range0.190592–5.152868eV. Ar9 RMS0.00760525eV. All narrow
regions remain in the metrics; the electron maximum0.94810eV includes displaced
piecewise-constant boundaries and is reported alongside the volume norm.
No source-time or thermal convergence follows. The enriched5ns cases below
will determine whether the measured mechanical correction survives a full run.

At the current continuation, poll the actual handles below; older live entries
elsewhere are historical. No result is complete until the full target is reached.

- **16195**: `reactive_material_argon_factored_heavy_5ns`, frozen
  `/private/tmp/frc-acoustic-factored-20261017`; complete20ns request with
  heavy-only enrichment and algebraically reordered sound/temperature/viscosity.
  The equivalent sound/temperature/Q expressions are now integrated into main.
- **61204**: `reactive_material_argon_scaled_heavy_5ns`, frozen
  `/private/tmp/frc-scaled-thermal-20261017`, same full20ns target with128-bit
  number/energy coordinate shift to retain the measured mixing heat.

- **31198 EXIT0**: `reactive_material_argon_e256_h32_scaled_10ns`, frozen
  `/private/tmp/frc-e256-scaled-20261017`; completed20ns with256 initial
  electron cells,32 initial heavy cells, variance return,128-bit scaling and
  no source-driven enrichment. Source SHA
  `64eb503a955d9ce992fd43a183cb62e44277ba9ec35ea86925eb8149921bde17`.
  Both intervals accept4869/3778steps; no source remains live in this root.
  The phase cap follows the declared maximum initial resolution. The completed
  comparison above with e128/h32 and the existing fine reference includes the measured
  arithmetic correction, not a pure single-parameter spatial-order study.

- **91194**: `reactive_material_argon_e128_full_uniform_h32_10ns`, frozen
  `/private/tmp/frc-full-uniform-electron-20261017`; corrected full20ns request
  with128 actual uniform initial electron velocity cells. Source SHA
  `1e2873c2ff590fbadb4d15c9ded1190e77579e63ee376587d78e68d0ab8e575a`.
  First10ns accepts2267steps,8469density cells,Erel−4.83442e-12; second interval
  is active. The original mixed-mesh9098 is terminal EXIT0, source
  `6bc02e83cb4ad514602f216aedf8c2027c5bc6304e7921f825155d91946d6370`.
  Fine-reference comparison selects new cases only when the full endpoint is
  successful. Do not substitute a prefix or the flawed setup for its correction.

Keep the old frozen roots intact. Compare complete 10/5 ns local temperature,
velocity, fronts and E/J when available. No 2.5 ns run yet: resolve the present
thermal/arithmetic errors before further mesh growth. The 20 ns planar component
is not the 1 m lab. Gas chemistry/collisions, heat, driven coupling, boundaries
and independent physical formation evidence all remain required.
