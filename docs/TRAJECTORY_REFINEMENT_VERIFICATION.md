# Refined trajectory rejection and independent density evidence

The full FRC objective and all 12 acceptance beads remain open. No unit tests
were added or run. This is a filled periodic contact verification fixture,
not a formation experiment or proposed laboratory operating point.

## Actual full-duration attempt

`outputs/trajectory_refined_t60_minus` requested the unchanged physical
`split_mode_evolution` protocol: F8/64,B32/64, epsilon=-1e-4,100us,60 Gauss steps,
absolute stage tolerance1e-12 and relative tolerance1e-10. The physical source is
`6a96344a64005476d0927d62ebaa772fb0c220dad584925b0363a661e3179114`.

It accepted 17 steps through28.333333us and rejected step18:
`Newton correction did not converge (GMRES info=5); no input modified`.
The terminal driver retains the accepted state, stage changes and iteration
histories. The exception text is not independent proof of in-memory atomicity.
Terminal raw-state comparisons are recorded separately in
`outputs/trajectory_refinement_diagnostics/primary_terminal_preservation.json`.

Session36840/PID89880 exited0 and was collected: the driver handles a numerical
failure by recording it. Exit0 therefore does not mean full-duration success.
The control explicitly records trajectory_failure=true and
full_duration_complete=false. Peak observed memory was4.439148GiB and minimum
free disk390.625MiB, within the8GiB/256MiB limits; no resource guard triggered.
This separates numerical rejection from resource interruption.

## Independent accepted-prefix density observation

New frozen `examples.trajectory_density_prefix` reconstructs every saved
endpoint through step17 at independent cylindrical quadrature16/20. Its source
is `a4f2c3743fe3441930c3c5be764dc236c45a2006d2db2b54b11cf4729aa66a71`.
The immutable input snapshot is
`outputs/trajectory_refinement_diagnostics/prefix17_manifest.json`.
Original summaries may later advance; all selected state/stage bytes are fixed.

It compares actual density against the independent radial1/N64 linear target
`n_base + epsilon*exp(lambda*t)*dn_reference`, lambda=8763.74925436694/s.
Initial representation error is retained. Saved coefficients are subtracted
before independent reconstruction to reduce background cancellation. This
observation constructs no coupled solver or induction matrices.

At28.333333us, order20 gives:

| Negative trajectory | Max nodal change / initial density | Density L2 error / growing linear signal | Actual density-change L2 / initial linear signal |
|---|---:|---:|---:|
| New F8/64,B32/64 split | 0.1602334652 | 790.2581025 | 1012.9945733 |
| Older F4/64,B16/64 guarded | 1.0392862e-5 | 5.6995976e-5 | 0.28181397 |
| Earlier F4/32,B8/32 split | 1.0374725e-5 | 3.8982815e-4 | 0.28179943 |

The new density error relative to the growing signal is1.51320e-5 at8.333us,
0.0147262 at16.667us,34.9739 at25us and790.258 at28.333us. The maximum nodal
density departure reaches16.0233%; independent order20 sampling reaches15.8652%.
Sampled minimum density is7.4908035e19/m^3. Positive density and accepted
nonlinear residuals do not establish accurate evolution.

Across these three prefixes, the maximum absolute order16/20 difference in
the reported relative L2 error is5.05108e-11. This is quadrature agreement for
that integral, not trajectory convergence. In the refined case, sampled
maximum fractional density change differs by8.59480e-4 between orders; maxima
are less resolved than the integral and do not share its agreement claim.

These are single-sign observations, not centered linear responses. The linear
reference is not an exact nonlinear solution at the large refined departure.
The older medium trajectory uses different induction arithmetic; source and
mesh differences remain explicit. No causal attribution to a particular
spatial mode, time integrator, arithmetic operation or physical instability is
established. No state was fitted, filtered, repaired or replaced.

Both full density observations completed with collected exit0. The separate
verifier confirms every non-runtime metric, sampled-field fingerprint and all
three raw time-series files match bitwise and have identical raw SHA256s:

- `outputs/trajectory_density_prefix17`
- `outputs/trajectory_density_prefix17_replay`
- `outputs/trajectory_density_prefix17_verification`

The primary observation's own cumulative peak was252.281MiB. Resource records
are final observations, not continuous guards. This small density observation
does not substitute for the full physical/conservation or trajectory replay.

## Replayed spatial localization

`examples.trajectory_density_localization` uses the same independently
reconstructed density fields and requires every original field/reference/weight
fingerprint to match. Its source is
`16d7c4b64a24f205bd439ca920b53de599f4e8dd06fe3b5194be022273ff0626`.
Both full observations have collected exit0, and the separate verifier proves
all non-runtime metrics and all three raw array files replay with identical
bytes. Evidence: `outputs/trajectory_density_localization17`, its `_replay`
and `_verification` directories. No additional physical trajectory was run.

The complete error is partitioned over fluid radial elements and axial cell
Fourier modes; no state or mode is changed. Relative radial-partition defects
are at most2.14092e-16 and discrete Parseval defects at most1.14701e-15.
These are identity checks, not a claim of physical convergence.

For the refined trajectory, low axial cell modes abs(k)<=4 carry99.7635% of
the squared error at13.333us,99.7891% at16.667us,99.7976% at25us and99.9051%
at28.333us. At the last endpoint, combined abs(k)=0,1,2 fractions are23.5262%,
43.6288%,24.6828%. Upper-quarter modes abs(k)>=24 carry only2.58258e-12 of
the final squared error. These are cell-translation modes with intra-cell
quadrature retained and possible alias folding, not an unaliased spectrum.
They do not support treating this trajectory departure as an axial near-grid
density error analogous to the separate high-radial-mode rate-noise study.

The final squared error is distributed over all eight radial elements:
6.6263%,16.1298%,21.6273%,22.0976%,17.8169%,10.8381%,4.2499%,0.6141% from
axis to wall. The largest volume contribution is the fourth element; this
differs from the largest nodal departure near the first internal interface.

At that first interface, r/R=1/8, the order20 density-jump L2 norm divided by
background trace density is2.27545%. The largest sampled interface jump over
the background reaches4.69502%; older medium/coarse prefixes reach only
4.83967e-10/3.89976e-10 as fractions. Jumps start at exactly zero. Surface
jump integrals use2*pi*r*dz and have different dimensions from volume norms;
their normalization remains separate. A trace jump is not itself a failed
discrete local conservation budget or proof of the instability's mechanism.

All values at both quadratures are retained. Tiny initial error fractions are
sensitive: the largest radial-fraction order16/20 change is0.0265357 at the
refined initial state, whose total relative density error is3.26e-9. For the
refined noninitial states it is at most3.55996e-6; from step8 onward it is at
most5.25928e-8. Across all prefixes the low-mode fraction changes by at most
1.20685e-5 and the sampled maximum jump fraction by2.37051e-7. Integral agreement
must not be transferred to initial fractional partitions or sampled maxima.

Primary localization peak memory was250.484MiB (final cumulative observation).
The scientific report `outputs/trajectory_rejection_report/verification.png`
was visually inspected: all four panels, axes, legends and limits are legible.
The PDF uses the same figure and was not separately inspected. Source, original
failed trajectory, independent observations and replays remain immutable.

## Active next action

The full original rejected protocol is now replaying in a fresh process, with
all original tolerances, initial data and stage/state recording:
`outputs/trajectory_rejection_campaign`. Parent unified session40737/PID93211,
child PID93212. Plan: `rejection_replay_plan.json` in the diagnostics directory.
This is a one-job plan; the earlier seven-job plan was not launched.

Poll the actual parent and collect its terminal exit; never restart on an
observation timeout. Confirm the same failed step, all17 state/stage files,
initial/final arrays and all non-runtime failure metrics before claiming a
reproduced rejection. A reproduced rejection is still a failed100us gate.
Then probe the original spatial RHS along the complete saved early error
direction, with amplitude controls and an independent-mode control, as specified
in TRAJECTORY_TANGENT_PROBE_DESIGN.md. This proposed probe is not implemented or
executed yet. Its purpose is to test whether the actual spatial derivative
amplifies the observed direction without using a time stepper. Opposite-sign,
same-grid arithmetic and time controls remain required. Do not lower bounds or
shorten the requested physical duration to obtain an apparent pass.

Full signed field/current/E, pressure, conservation and perturbation-energy
audits remain unfinished, along with local budgets, heated exact vacuum,
transport/ionization, real prescribed coil/end boundaries, converged formation
and the experimental UI. No3-D stability or particle-confinement claim follows.
