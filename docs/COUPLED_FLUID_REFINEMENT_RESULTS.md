# Full fluid-resolution comparison — 14 October 2026

The original full coupled RHS now has completed F8/F16/F32 physical comparisons
at fixed B32/64. Both fine probes, both physical audits and the cross-grid
comparison have complete fresh replays. The original equations, physical
masses, wall constraints and induction tolerance 1e-10 were unchanged.
This is a completed instantaneous-rate study; it has not cured or advanced the
failed F8 trajectory.

## What the result changes

For the saved error direction at amplitude 1e-5, q20:

| Measured quantity | F8 | F16 | F32 |
|---|---:|---:|---:|
| Physical density projection, /s | 950294.3 | 952121.7 | 955179.6 |
| Physical radial-velocity projection, /s | 951737.1 | 655903.8 | 504650.2 |
| Physical Btheta projection, /s | 950317.7821 | 950317.7824 | 950317.7818 |
| Physical net power, positive sign, W | -4.53698e-10 | 2.52026e-11 | 1.05341e-12 |
| Absolute net power / sum of absolute component powers | 0.0963472 | 0.00595800 | 0.000247617 |

These projections are directional measurements, not eigenvalues. The F16/F32
radial-velocity responses have orthogonal fractions 0.998922/0.999774 relative
to the saved direction, so its scalar projection especially cannot describe
the response. Native grid net powers remain near roundoff while the separately
integrated physical power exposes an error that decreases with refinement.
No energy correction was applied. This does not close local or trajectory
energy balances.

Actual F16→F32 physical rate differences divided by F16 rate norms are 0.830318
for density, 1.21560 for radial velocity, 0.0901282 for electron entropy,
1.49970 for electron pressure, 2.376e-8 for Btheta and 1.190e-8 for Jz.
The strong fluid response is not converged; the rapid magnetic response is
almost unchanged by this particular refinement. The fixed saved direction can
contain interface jumps, so a larger strong derivative is not by itself a
diagnosis of the full trajectory's failure.

The independent smooth Fourier/Bessel mode remains a separate control. F32
rate errors include density 4.27815e-4, radial velocity 4.37484e-5, Br 0.0139750,
Jr 0.0402535 and Jz 0.00508541. Radial fluid refinement alone does not resolve
the remaining magnetic/current error. Finite-amplitude and base-roundoff
effects remain in the full saved metrics; no acceptance is inferred from the
selected 1e-5 row.

## Evidence and next simulation

Raw data and full metrics: `outputs/coupled_fluid_probe_F16`, `_F32`,
`outputs/coupled_fluid_audit_F16`, `_F32`, and
`outputs/coupled_refinement_comparison`, each with a separate replay and
verification directory. The cross-grid comparison retains all 17 physical
components, three amplitudes, both directions and quadrature16/20. Its source
identity is `905f775324e359d8b57d66664f6db654c7fd3fdf8f0c41f5e4d315d2d528c8f4`.
The F32 audit pair and comparison pair completed with collected exit codes0.

Following the user's simulation-first instruction, the next action is an
actual full 100 µs evolution at F32/64,B32/64 with the existing analytic initial
mode, epsilon=-1e-4 and60 steps. It keeps the original1e-12 absolute and1e-10
relative stage checks and8GiB/256MiB resource bounds. No new verification
framework, artificial diffusion, floor, state repair or modified mass is added.
Compare its saved evolution directly against the independent baseline before
choosing a solver change. Full acceptance still requires physical accuracy,
appropriate signed/time controls and the remaining project gates.

## Full evolution outcome and time-step refinement

The F32/64,B32/64 60-step run is terminal: it rejected the first step after
816.244 s, with a Gauss line-search residual of 1.8417582888632157e-12.
No time step was accepted. Peak own RSS was 5,282,955,264 bytes; no resource
limit triggered. The campaign process returned exit 0 through its recorded
numerical-failure path; this is not a successful physical trajectory.
Evidence: outputs/trajectory_fluid32_t60_minus and its terminal control in
outputs/coupled_refinement_diagnostics/fluid32_t60_minus_control.

A separate full 100 us run with 120 steps halves dt to 0.833333 us, keeping
all physical inputs and residual bounds. It runs under the existing driver
at outputs/trajectory_fluid32_t120_minus. This is a time-step control, not
an accepted cure. Compare actual saved evolution with the independent mode
when accepted states are available. No failure replay, new wrapper or unit
test was added. See GOAL_STATUS.md for the current worker.

The 120-step control also rejected step1, at residual1.0288129410491659e-12
after531.094s, peak5,467,111,424 bytes. Actual session55929 exit0 was collected;
again zero accepted steps and no resource stop. Halving dt nearly halves the
stalled residual, consistent with a rate-evaluation accuracy limit but not
sufficient to identify its cause. The next actual full100us trajectory uses
240 steps (dt0.416667us), all other inputs and bounds unchanged.

## Independent shorter-wavelength comparison

Using the existing original-equation Fourier/Bessel operator, radial modes
32,64,128 have maximum real eigenvalues12513.4633,13325.6316,23551.9380 /s
at128 axial points. Doubling axial points64→128 changes these by1.406e-7,
2.730e-8,1.735e-11 relative. Leading eigenvalues are real; relative divergence
residuals are1.762e-9,1.745e-9,1.059e-9. This extends the mathematical baseline
without changing production physics. Raw operators, eigenvalues and vectors
are in outputs/contact_short_radial_spectrum; source closure is included.
Actual process5184 exited0; peak ownRSS243,351,552 bytes. No replay or unit tests.

These sampled continuum growth rates remain far below the approximately
950,000 /s production saved-error projection. A finite set of radial modes
is not a bound on all modes, and a directional projection is not an eigenvalue;
this does not establish a numerical cause or certify a physical instability.
The next useful evidence remains the finer fluid trajectory itself.
