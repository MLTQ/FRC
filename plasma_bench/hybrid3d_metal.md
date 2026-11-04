# hybrid3d_metal.py

Optional MLX 0.32.4 custom Metal kernels for the actual hybrid device solver.
Import is lazy; CPU runs do not require MLX. An explicitly requested GPU
raises if unavailable; there is no silent fallback to CPU.

- `MetalKernels.transport` fuses upwind internal-energy face fluxes and wall
  outflow rate, using float32 arithmetic and immutable staircase geometry.
  It synchronizes and returns float64 host arrays. Conservative divergence,
  compatible pressure work and stored electron energy remain CPU float64.
  The same kernel returns outgoing transport rates for the independent
  positivity CFL; a 1e-6 conservative margin covers its float32 arithmetic.
- `MetalKernels.moments` is an experimental normalized float32 atomic CIC
  deposit. Coordinate normalization now runs on GPU; density and velocity scales
  keep atomic values near unity. No mass, temperature or energy is fitted.
  Atomics may be nondeterministic and do not have CPU64 precision.
  Stored particle positions remain CPU64; gather coordinates are cast to
  float32. `prepare_moments`/`moments_device` expose a separate resident-input
  timing; it excludes upload/host collection and is not a whole-step speedup.
  Equal destination keys in groups of eight SIMD lanes are combined before
  atomic writes. Other keys use the general atomic path; no particle order,
  geometric symmetry or spatial sector is required for correctness. Speed
  depends on particle locality, so startup timing cannot imply late-shot speed.
- Timing counters include casts, uploads, execution, synchronization and host
  conversion. JIT warmup must be separated in measurements.

This is a mixed-precision acceleration prototype, not a validated all-GPU
solver. CPU64 particle positions/velocities, Faraday curl, exact vacuum solve
and circuit are retained. Matched physical replays and independent analytic
comparisons are required; microkernel speed is not a full-shot speedup.
GPU execution may require a non-sandboxed process on macOS.

Measured 1.29-million-ion startup deposition after the transport-CFL fix:
CPU two-thread 21.59 ms, Metal including transfers 6.79 ms (3.18×),
resident-input kernel/sync 1.80 ms (12.00×). Complete startup throughput
improves 1.149×; a coupled 0.75 µs coarse replay is essentially equal in
runtime. These are different ratios. See `docs/HYBRID_3D.md` for evidence
and physical comparisons; no late-shot performance claim follows.
