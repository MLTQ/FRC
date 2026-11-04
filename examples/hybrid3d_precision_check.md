# hybrid3d_precision_check.py

## Purpose
Bounded full-device physical replays isolating the two Metal float32 blocks,
for frc-ibl.11.10.1. It reuses existing initialization, physics, diagnostics and
analytic CIC comparison. No production solver edit, unit suite or new campaign.

## Components and contracts
- `main`: new output only, freezes both Python packages and hashes through the
  existing campaign helper, then executes the worker from that immutable copy.
  Input formation data remain an explicit absolute path outside the source copy.
- `worker`: sequential CPU64, identical CPU64 repeat, transport-only,
  deposition-only and combined runs. Default actual 0–0.1 µs, 10 mm isotropic,
  0.5 ns, eight ions/cell, 5% cutoff, same seed-0 initial fill and coupled banks.
  At most one extra physical replay runs; current long jobs are never changed.
- Deposition-only installs `MetalKernels` into the existing `_metal_moments`
  extension point after ordinary CPU initialization. Electron transport stays
  CPU64. All cases begin from identical CPU64 initial density/fields/particles.
- `snapshot`: native float64 fields, pressure/energy, full particles and flux
  arrays at initial/middle/end. No viewer float32 truncation masks a difference.
- `moment_probe`: actual native initial/endpoint particles deposited repeatedly
  on Metal, compared with CPU64 and with input-rounded CPU64 accumulation.
  Exact integrated particle number/flux provide an independent global reference.
- `compare`: matched times against CPU; dense/thin/empty support uses mean
  density and actual gas dual volumes. Unweighted node norms are not energy
  integrals; CPU is a reference, not physical truth. Particle identity is not
  paired if absorption changes the population count.

## Equations and limits
Unchanged kinetic-ion/quasi-neutral-electron hybrid with Faraday/RK3, physical
Spitzer resistivity, declared density cutoff and explicit cutoff/empty-node
heat sinks. Exact geometric vacuum, absorbing walls and source work are retained.
No mode/symmetry filter or energy/temperature adjustment. Electron internal-energy
upwind transport is first order/diffusive, as documented by its module. This
comparison isolates arithmetic, not the physical validity of those loss rules.
The short initial-fill window does not establish late-time statistical or 3D
convergence. Metal atomics may differ between identical calls.

## Use
`NUMBA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m
examples.hybrid3d_precision_check outputs/NEW --input outputs/fine_uniform_h2`.
GPU-capable execution is required. `progress.json` and `worker.log` report actual
progress/failure. A failure is retained rather than rescaling or silently retrying.
Saved diagnostics use the existing campaign's max(5 J, 1% of source work)
accounting guard; exceeding it stops that owned comparison without state repair.

Optional `--modes split_deposition split_combined` adds only the new CPU64-cell /
Metal32-fraction CIC prototype. `--reference EXISTING_CPU_RUN` compares against
the existing native CPU replay without rerunning it. Native times and grid shapes
must match; initial fields/particles are compared as well. Each addition gets a
new source snapshot/output; an already running comparison retains its old code.
The reference can also be an existing source-frozen long CPU run. Its actual
matched `snap_*.npz` times are read, full fields compared and precision labelled
`saved_float32`; no incomplete ion display sample is paired as full particles.
This avoids duplicating the long CPU simulation, but small differences are
limited by its saved precision and cannot prove native float64 agreement.

## Measured evidence (2026-10-22)
`outputs/hybrid3d_precision_isolation_fine_v1/`: all five cases completed actual
0–0.1 µs. CPU repeat is bitwise identical, including signed zero/dtypes, in every
saved native field/particle/flux array. Its CPU baseline also exactly reproduces
the long control's float32 B/n/pe at 0.1 µs. Response-field relative differences
are 7.16e-11 transport-only, 4.869e-6 deposition-only and 4.869e-6 combined.
Density/pressure differences with original deposition are 1.897e-6/1.957e-5.
The pressure-work ledger residual is about -0.003117 J against 24.598 J source
work in every case; numerical heat removal is 0.021813 J. Those small early
budgets do not validate the later large numerical cooling or identify which
late trajectory is physical. Region norms are unweighted numerical diagnostics.

`outputs/hybrid3d_precision_split_fine_v2/` completes both revised deposition
cases over the same actual interval. Combined response difference falls from
4.86900e-6 to 5.61287e-9 (867.47x); density is 1.08054e-7 and pressure
9.84041e-9. All 21 native states across the seven cases are finite and their
source snapshots remain intact. The source-frozen 0–1 µs combined follow-up is
`outputs/hybrid3d_precision_split_to1us_v3/`, with the existing long CPU archive
as reference. Exercising that archive adapter on the completed 0.1 µs prototype
gives response difference 2.19966e-8 instead of native 5.61287e-9: this measures
the limitation from saved-float32 reference fields, not additional evolution.

That 1 µs extension has completed all 2,000 steps with finite native states
and intact source hashes. At 0.5/1 µs, revised response differences are
3.21292e-5/0.0141901, against original GPU 0.00196788/0.0950132: reductions
61.25/6.70x, so the short-window 870x factor does not persist. Dense-region
field difference at 1 µs is 1.08908e-4; low-density differences still grow.
Source work 1,001.529 J, ledger -0.26694 J and numerical removal 157.169 J
show the cooling is not cured. `summary.json` retains actual frame paths,
matched budgets and comparison precision. Next is a same-seed corrected repeat
before any longer adoption; production defaults and live jobs are untouched.
