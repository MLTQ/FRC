# hybrid3d_gpu_check.py

Runs the actual optional Metal backend, independent physical comparisons and
synchronized whole-device timings. It is not a unit-test suite.

- `independent` uses the existing moving uniform-pressure contact and
  constant/Spitzer force-free decay against analytic translation and an
  independent continuum ODE, at 16/32/64 cells. These are periodic fixtures,
  not formation or wall-loss experiments.
  Transport-only and transport-plus-CIC are recorded separately.
  An independent uniform-density/current deposition fixture also checks
  ordered/shuffled particles and an incomplete final SIMD group against
  exact number/current integrals, not only agreement with the CPU kernel.
  Gaussian transport at outer CFL≈1.5 is checked separately at32/64/128cells:
  the exact solution is periodic translation; the operator must subcycle,
  remain positive and conserve energy while its diffusion decreases.
- `device` advances the same true initial-fill coupled shot with CPU64,
  Metal transport and Metal transport plus CIC. It excludes the first ten
  steps/JIT warmup from timing. All conversions, GPU synchronization, CPU
  fields/vacuum/circuit and particle work remain in the measured step time.
  Endpoint snapshots and budgets are retained. The default 40-step interval
  is startup, not a prediction of full-shot or formed-plasma speedup.
  `--modes 1 2` can add only the GPU replays after a saved CPU measurement;
  `device.json` retains earlier rows rather than silently replacing them.
- `particle_block` times real full-resolution device particles through CPU
  and normalized GPU CIC, reporting node errors and independent deposited
  particle-number closure. It is a kernel timing, separately labeled.
  The resident-input timing excludes uploads and CPU collection; it cannot
  be reported as the speedup of the current complete mixed-precision step.

Use `NUMBA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 .venv/bin/python -m
examples.hybrid3d_gpu_check outputs/hybrid3d_gpu_check`. macOS sandboxed
processes may not see Metal; use authorized GPU-capable execution. Evidence
is ignored output JSON/NPZ. MLX 0.32.4 and mlx-metal 0.32.4 are optional
dependencies; CPU production remains the default and active source-frozen
campaigns are not switched by running this example.
The entry point records source hashes; use fresh output directories for
changed kernels so evidence is not mixed across revisions.
