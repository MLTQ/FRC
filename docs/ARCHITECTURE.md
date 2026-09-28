# Architecture and stack decision

Python 3.11+, float64 NumPy arrays, SciPy elliptic integrals and sparse LU,
Matplotlib/ContourPy visualization and topology, pytest verification.
The initial environment is macOS ARM64 / Python 3.13, not the eventual RTX 4090
machine. CPU correctness comes first; no accelerator is assumed or required.

NumPy/SciPy provide the required special functions and a sparse direct solver
without writing numerical infrastructure. Geometry-dependent loop bases and
diffusion factorizations are reused. For this modest axisymmetric problem, a GPU
transfer/JIT layer would add complexity before we have a performance bottleneck.
JAX can support future differentiable/JIT array kernels, Numba can accelerate
profiled local loops, and PETSc becomes attractive for larger implicit solves.
None is required by a public interface. Revisit after measured profiling and the
Stage 2 solver choice, on the actual target GPU operating system.

References: [SciPy cached factorization](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.factorized.html),
[JAX supported installation platforms](https://docs.jax.dev/en/latest/installation.html).

```
geometry + typed configuration
      |         waveforms -> coil definitions
      |                          |
      +---- analytic vacuum field basis
      |                          |
      +---- plasma solver protocol (Stage 0 / Stage 1 / future MHD)
                                 |
                     result arrays + metadata
                                 |
                    diagnostics / topology
                                 |
             plots / time slider / sweeps / saved NPZ
```

Modules have one responsibility and adjacent .md contracts. Arrays use `(nr,nz)`
and histories `(nt,nr,nz)`. SI fields and axis conventions are shared; coils do
not import plasma solvers. Frozen dataclasses validate geometry, material and
source configuration. Waveform objects supply value and derivative; callbacks
are script-only and rejected by reproducible serialization.

Results include full configuration, source/version fingerprint, dependency
versions, actual time coordinates, warnings, units, boundary conditions, model
stage and numerical method. Saved archives are numeric NPZ without pickle;
metadata is JSON. JSON configurations and a CLI keep runs reproducible without
a notebook. Sweeps replace typed configurations, not hidden globals.

The developmental dynamic kernel is `mhd.py`, with explicit modules for annular
mesh geometry, conservative state, Rusanov fluxes, compatible induction and
analytic benchmarks. It evolves density, poloidal momentum, total energy and
corner flux with SSP RK2. Its first-order spatial dissipation is measured.
It fills a verification cylinder with fluid, with a conducting radial wall and
periodic z. The original Stage 2 formation entry point still refuses execution
until moving-plasma/vacuum and source-work coupling are verified; it does not
silently substitute these verification boundaries for the physical experiment.
