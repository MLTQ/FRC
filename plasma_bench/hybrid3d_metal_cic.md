# hybrid3d_metal_cic.py

## Purpose
Isolated accuracy prototype for Metal CIC, owned by the numerical lane
frc-ibl.11.10.1. It does not replace the existing kernel or alter running shots.

## Components and contracts
`SplitCoordinateMetal(MetalKernels)` retains the existing density/velocity unit
normalization, float32 weights/velocities, eight-lane combine and atomic sums.
Only position handling changes: evaluate s=(x-lo)/d and floor(s) in CPU64,
send integer cell indices and float32 fractions to the GPU. Whole-cell indices
and subtraction of world-coordinate origins no longer consume fraction bits.
This is the same mathematical CIC shape used by the CPU64 gather, with explicit
remaining fraction/weight/velocity/accumulation rounding; no weight, charge,
energy, temperature or mode is fitted. Periodic wrapping is unchanged.

`prepare_moments` and `moments_device` preserve the base API; the inherited
`moments` includes host coordinate arithmetic, uploads and synchronization in
its measured time. It sends an additional three int32 values per particle.
This may reduce speed; actual whole-step comparisons are required.
`split_coordinates` fuses CPU64 normalization/floor/fraction into independent
Numba particle writes, without fast-math, avoiding full-sized float64 temporary
coordinate arrays. It has no reduction or particle-order assumption.

## Verification and scope
Use independent uniform density/current/number, ordered and shuffled/partial
groups, manufactured linear-electric-field power and actual coupled startup
replays. Compare full native fields, particle/energy budgets and transport-only
and CIC-only cases. Float32 atomics remain nonrepeatable. A smaller early error
does not establish long-time stability or physical vacuum-edge losses. All
production imports keep the original implementation until evidence supports
an explicit opt-in; current source-frozen campaigns are never changed.

## Measured evidence (2026-10-22)
The independent uniform density/current, shuffled and partial-group comparisons
pass. On 1.29M native evolved particles at 0.1 µs, raw density/flux differences
from CPU64 fall 1.68e-6/7.46e-6 to 1.11e-7/1.39e-7. Repeated GPU atomics still
differ at about 8e-8. The independent quadrupole-work error remains about
2.4e-9 of absolute particle power; it is not improved by this coordinate change.
Evidence: `outputs/hybrid3d_cic_coordinates_endpoint_v2/`.

The actual 10 mm, 0.5 ns, 0–0.1 µs deposition-only coupled replay gives response
field difference 5.60e-9 versus 4.87e-6 for original CIC (~869x smaller), density
1.08e-7 and pressure 9.80e-9, with unchanged physical equations and initial state.
The combined transport/deposition comparison completes with 5.61e-9 response
difference, 1.08e-7 density and 9.84e-9 pressure in
`outputs/hybrid3d_precision_split_fine_v2/`. A 0–1 µs extension uses the existing
CPU control's actual fields; its float32 save precision is disclosed. It has
now completed: response difference is 3.21e-5 at 0.5 µs and 0.01419 at 1 µs,
61.25/6.70x below original deposition. Numerical removal remains 157 J; the
initial 870x reduction is not a long-time factor. Next is a same-seed corrected
repeat to measure residual atomic sensitivity before longer adoption. Long-time
sensitivity remains open; the existing live campaigns retain the old kernel.
