# Shared bulk transport and physical thermal verification

This checkpoint advances actual conservative density and both species log
entropies in prescribed moving flows. **370 tests pass in 126.72 s; 14 rate
cases and 10 complete trajectories are independently audited.** Six diagnostic
pressure-force probes have bitwise array replay. Another 24 linear controls
isolate and verify the correction of an initially unstable face reconstruction.
The full momentum/current/magnetic evolution has not adopted these components.

Equations, state, boundaries and numerical limitations are in
[BULK_THERMAL_NUMERICS](BULK_THERMAL_NUMERICS.md). The cylinder has a regular
axis, impermeable radial wall and periodic axial coordinate. This smooth,
filled-cylinder verification is not a formation or free-edge experiment.
There are no added viscosity, resistivity, density floors or energy resets.

## Independent continuum rates

The characteristic reference uses R=0.035 m, axial half-length 0.12 m,
n0=1e18 m^-3 and speed 2000 m/s. Electron/ion gamma are 5/3 and 1.4. Its
continuity and pressure equations are checked independently by complex-step
derivatives. Physical norms use separate root-product reconstruction; the
production moment/primitive implementation is not its own reference.

Degree-five relative volume L2 errors, with nz=2 nr:

| nr | Density rate | Electron log-entropy rate | Electron internal-energy rate |
|---:|---:|---:|---:|
| 8 | 1.35930e-5 | 6.80113e-3 | 1.19946e-3 |
| 16 | 2.08821e-7 | 1.32368e-4 | 2.28213e-5 |
| 32 | 3.24703e-9 | 2.24888e-6 | 3.86994e-7 |
| 64 | 5.06720e-11 | 3.58944e-8 | 6.17466e-9 |

The finest first-annulus electron entropy-rate error is 3.19378e-8.
Degree-three, physical quadrature, density 1e16/1e20 m^-3, reverse flow,
constant entropy and stronger density contrast controls are retained. Constant
entropy has exactly zero rate. Stronger contrast gives electron entropy-rate
error 1.31003e-4 at nr32; this is not silently grouped with the weaker case.

Across all 14 rates, maximum normalized pressure-work backward error is
2.34e-15, scalar solve residual 1.95e-13, independent local-energy agreement
4.37e-15, physical entropy-rate defect 1.78e-15 and weak continuity defect
8.78e-16. These are discrete identities and replay checks, distinct from the
continuum accuracy in the table.

## Actual evolving density and two species temperatures

Uniform translation runs a complete axial transit of 120 microseconds.
Compression uses v_z=V sin(kz), integrated to |k V t|=ln(2), or
13.238136009 microseconds. Each RK4 stage computes particle momentum from its
current numerical density and prescribed velocity; exact evolving density is
never injected. Five physical snapshots and every accepted-step diagnostic
are saved. The reference satisfies the prescribed-flow thermal equations,
not the unforced momentum equations.

Final relative physical volume L2 errors:

| Flow | nr | Density | Electron temperature | Ion temperature |
|---|---:|---:|---:|---:|
| Compression | 8 | 1.60995e-3 | 1.00875e-3 | 5.06501e-4 |
| Compression | 16 | 7.05037e-5 | 4.28221e-5 | 1.92981e-5 |
| Compression | 32 | 1.37162e-6 | 8.26524e-7 | 3.73944e-7 |
| Full translation | 8 | 2.81597e-5 | 2.19714e-5 | 1.16304e-5 |
| Full translation | 16 | 4.66067e-7 | 3.62877e-7 | 2.04565e-7 |
| Full translation | 32 | 8.33719e-9 | 6.57306e-9 | 4.39992e-9 |

Initial representation errors remain in the saved norms. For nr32 translation,
initial density and electron-temperature errors are 9.82e-11 and 6.78e-10;
the final table is not a claim of initially exact representation. Entropy-wave
norms remove the background, rather than hiding wave error under a large mean.
Final electron/ion wave-relative errors there are 1.13e-8 and 9.43e-9.

At nr16, halving the compression time step changes final density by 7.55e-10
and electron temperature by 5.83e-10 in independent physical relative L2.
This separates spatial error from the tested time error. Quadrature order 16,
reversed compression and degree-three controls complete the ten-case protocol.

Particle-number drift stays below 2.23e-16. Physical entropy is bilinear in
evolving density and entropy coefficients, so RK4 need not preserve it exactly:
maximum drift is 2.59e-9 in coarse translation, falling to 2.54e-12 for the
finest ion translation. Electron entropy can be much closer to roundoff in this
particular symmetric fixture; that is not an all-state conservation theorem.

The physical internal-energy change is compared with an independently integrated
pressure-work ledger. No temperature is reset from that ledger. Its maximum
discrepancy divided by initial species energy is 1.32e-8 across the protocol.
For nr32 compression it is 4.20e-13 (electron) and 4.09e-13 (ion). At nr16,
time halving reduces these defects by factors 17.41 and 15.25. Translation has
nearly zero continuum net work: a small initial-energy-normalized defect must
not be mistaken for a small error relative to that tiny energy change.

## Mechanical coupling is still a separate requirement

The six pressure probes compare a candidate conservative stress with independent
continuum pressure gradients and the actual reconstructed kinetic-energy
derivative. This candidate is diagnostic; it is not a momentum update.

| nr | Pressure-force relative L2 | Mechanical/thermal work backward defect |
|---:|---:|---:|
| 8 | 4.29435e-4 | 4.56114e-7 |
| 16 | 9.67024e-6 | 2.18128e-8 |
| 32 | 1.67948e-7 | 3.92847e-10 |
| 64 | 2.69703e-9 | 6.48018e-12 |

The second norm divides by uncancelled power, not by the nearly cancelled net
power. The finest signed discrepancy is -1.10429e-9 W. Reverse flow flips its
sign. Uniform analytic pressure still produces a false force of 2.57208e-5
Pa/m RMS at nr32, including 1.61502e-6 Pa/m in the first annulus. Thus the
candidate converges accurately but is not an exact pressure-work adjoint or
exact uniform-pressure balance. Simply treating the thermal face-work
coefficient as an acceleration would not solve conservative momentum balance.

## Failure preservation, storage recovery and audit scope

The [initial attempt](BULK_THERMAL_INITIAL_ATTEMPT.md) is preserved under
`outputs/*_initial_attempt`. An inappropriately biased normal-face fit caused
positive advection growth and nonpositive density in full translation. A
centered integer-index face fit corrects the cause; neither density floors nor
added diffusion hide it. The 24 independent face controls include both signs,
degrees and mesh sizes. Centered growth is below 4e-14 per transit at measured
precision. Original nr32 translation was deliberately stopped after diagnosis;
that interrupted run is not counted as completed.

The Mac subsequently ran out of disk during corrected result serialization.
Removing regenerable project caches restored local writes. The original worker
finished with a storage error; only absent cases 2, 4, 6, 7 and 8 were explicitly
recomputed, with unchanged scientific source. Already saved cases were retained.
The original log, separate recovery helper and `disk_recovery.json` document
this event. It is not a numerical retry or a repaired scientific state.

All 14 rate arrays were replayed bitwise. For all ten transports the audit
replays initialization and the last accepted RK4 step bitwise, checks the full
positive/finite accepted history and work increment, and independently samples
all five snapshots at quadrature order 18. It does not claim bitwise replay of
every earlier step. Independent snapshot integral agreement is at most 1.67e-15.
All six pressure-probe arrays also replay bitwise. The report checks raw and
audit identities, and its final six-panel figure was visually inspected.

Calculation source: `a74e877fe9dd2030a6957d717d87dc1a1e96c4866fdb92ba29e5b76c157fcd2b`.
Initial failed source: `a90606b784ad3ec01fc3eda3897b60330bac775c8f217fbe6bf650e69100fbc5`.
The previous log-thermal, variable-current, compatible-wave and smooth-moving
calculation identities are unchanged. Source archives, arrays, audits, logs,
report and checkpoint are under `outputs/bulk_thermal_*`; diagnostic mechanics
and face evidence are under `outputs/bulk_pressure_work_probe` and
`outputs/particle_face_diagnosis`.

![Measured bulk transport convergence](../outputs/bulk_thermal_verification/verification.png)

Next derive and verify the full common moving state: conservative momenta and
angular momentum, true physical kinetic-energy gradients, velocity/vorticity,
pressure/Maxwell/drift stresses, current and enthalpy/Poynting work. Then test
independent parent-species rates and complete nonlinear wave/contact evolution.
The old nonlinear solver's known errors remain; this checkpoint does not replace
them. Heated-edge/exact-vacuum coupling, physical transport/ionization, realistic
ends/coils, formation convergence and an accurate interactive experiment remain
unfinished. The overarching FRC objective remains open.
