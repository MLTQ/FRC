# Completed full rejection replay and spatial tangent — 2026-10-12

The subsequent complete term decomposition and one-sided current-flux audit
are recorded in SPATIAL_PARTS_VERIFICATION.md. The next-action paragraph below
records this tangent checkpoint's original handoff and is superseded there.

The negative F8/64, B32/64 trajectory's entire original protocol now replays:
17 accepted steps through 28.333333 µs, followed by step-18 Newton GMRES info=5
rejection. All 19 raw initial/final/state/stage files and non-runtime metrics
match exactly. The campaign and both physical workers are terminal, exits
collected. Replaying a rejection still fails the requested 100 µs gate.
See `trajectory_refined_t60_minus_replay_verification` and
`trajectory_tangent_diagnostics/rejection_completion.json`.

The full original spatial RHS was then applied to both signs of the complete
step-10 error direction and the independent continuum tangent, at normalized
amplitudes 1e-4, 1e-5 and 1e-6. No time integrator, filtered direction, changed
mass, floor, diffusion or relaxed bound was introduced. The declared fixed
coordinate metric has 22208 components. The nonzero base coordinate-rate
maximum is 1.1573374e-7/s and is retained explicitly.

All six signed cases passed the original 1e-10 induction bound; maximum signed
residual was 5.91859e-12. Actual passed state/field fingerprints stayed unchanged.
All 26 raw files and complete metrics reproduce in a fresh process. The
independent physical audit at quadrature 16/20 also completely replays, with
one identical raw integral file. All four probe/audit workers are terminal,
exits 0 collected, resource guards passed. Maximum audit memory was
2681962496 bytes; minimum free disk across these jobs was 366592000 bytes.

| Direction | Normalized amplitude | Coordinate projection (/s) | Orthogonal fraction |
|---|---:|---:|---:|
| Saved error | 1e-4 | 950740.878742 | 0.003730800 |
| Saved error | 1e-5 | 950740.874122 | 0.003730776 |
| Saved error | 1e-6 | 950740.870052 | 0.003730776 |
| Independent mode | 1e-4 | 8763.717808 | 0.000415950 |
| Independent mode | 1e-5 | 8763.717810 | 0.000415950 |
| Independent mode | 1e-6 | 8763.717899 | 0.000415971 |

These are directional projections, not certified eigenvalues. At q20 and
h=1e-5, independently reconstructed saved-error density projects at
950294.320239/s with orthogonal fraction 0.004416805. Electron/ion pressure
projections are 950229.990109/s and 949461.786879/s. A rapid direction therefore
exists in the actual spatial RHS without time integration, but its physical
versus numerical origin remains unresolved.

The physical-RHS derivative and conservative-Jacobian pushforward remain
distinct because the base has a small stationarity defect and reconstruction
is nonlinear. Their largest saved-error component difference decreases from
9.32686e-8 to 9.32364e-12 as h decreases. E is an algebraic response, not E_t.
Maximum relative q16/20 projection difference is 1.42018e-14; this scope does
not include pointwise maxima. See `quadrature_comparison.json`.

Independent-mode controls retain unresolved accuracy limits. At h=1e-5,
density rate error is 0.000427815, Br error 0.01397495, radial-current error
0.04025354 and azimuthal-current error 0.01110179. Adjacent-amplitude radial
current differences reach 0.000875113; no amplitude-convergence claim follows.
Saved-error amplitude differences reach 3.33422e-6 in radial current.

Next: decompose the unchanged spatial RHS along this saved error into physical
bulk, current, pressure/species, induction and wall-reaction contributions;
reconstruct the complete original sum and compare independent continuum and
interface terms. Do not remove a term as a proposed fix without proving a
defect. Opposite-sign full trajectories, same-grid induction-arithmetic
controls, time refinement and complete physical/local conservation remain.
The old seven-job plan was never launched and must not duplicate this replay.
All 12 overarching acceptance beads remain open. No unit tests.
