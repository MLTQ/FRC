# Established plasma engines: bounded feasibility comparison

Status: bounded comparison complete, 2026-09-29. Production physics and the Stage 2
formation gate are unchanged. The user authorized an isolated evaluation on
`m@192.168.0.202`, while preferring to retain existing work if migration does not
save substantial effort. This comparison does not select different physical
assumptions merely to obtain an attractive plasma animation.

## What is already reusable

| Existing capability | Retain if an external engine is adopted |
| --- | --- |
| Typed chamber, plasma, coil and waveform definitions | Experiment specification and an engine adapter |
| Analytic circular-loop flux and induced electric fields | Source construction and independent field checks |
| Topology, probe, scale and validity diagnostics | Common postprocessing, with explicit conversion conventions |
| Analytic wave, force, diffusion and equilibrium references | Independent verification suite for either engine |
| Exact massless magnetic exterior and source-work accounting | Current quasistatic reference; different equations from full Maxwell vacuum |
| Bounded MHD, Hall, electron/ion thermal and finite-mass components | Research references, retaining their documented limits |

The complete moving, heated, finite-mass plasma with a physical free edge,
real ends and spatial coil drive is **not implemented as one verified model**.
The current count of passing component tests cannot be converted into a credible
percentage completion or calendar estimate for that integration.

## Candidates and concrete compatibility findings

### WarpX

Released source: `26.09`, commit
`0c62c75e53a9ad08241535444bd7e53fd1deba88`.
AMReX `a52ca73324ac2c7b65ec04f131e6df99eec9c576` and pyAMReX
`9c90b065078badca6900ead74ee3ca482acfe186` are the fetched build dependencies.
The isolated build enables 1D, CUDA architecture 89, Python and double precision;
MPI, embedded boundaries, FFT, openPMD and QED are disabled for these fixtures.
It is not an RZ or 3D performance measurement.

Directly useful capabilities include kinetic ion trajectories, full kinetic
electromagnetics, GPU execution, a Python interface, and prescribed vector
potentials for changing external fields. The source already includes an RZ
plasma-cylinder compression example. Actual FRC stability use is reported by
[TAE/LBNL researchers](https://meetings-archive.aps.org/dpp/2023/tp11/105/).

The [released cylindrical example](https://github.com/BLAST-WarpX/warpx/blob/26.09/Examples/Tests/ohm_solver_cylinder_compression/inputs_test_rz_ohm_solver_cylinder_compression_picmi.py)
uses an initially pressure-balanced cylinder, periodic axial geometry, an
electron-pressure closure, density cutoff, prescribed resistance and additional
hyper-resistivity. It is not our initially preionized argon pulse-formation
problem, and a short regression run of it would not prove formation accuracy.

The [hybrid electric solve](https://github.com/BLAST-WarpX/warpx/blob/26.09/Source/FieldSolver/FiniteDifferenceSolver/HybridPICSolveE.cpp)
divides its Hall/pressure electric force by `max(rho,rho_floor)`. Its optional
Holmstrom branch removes those terms in low-density cells and retains resistive
magnetic evolution. That regularization is not our algebraic, current-free
magnetic exterior. A denominator floor does not itself load massive particles,
but it still changes the field equations; these distinctions must remain explicit.
Hybrid electrons have no inertia. Full electromagnetic PIC is a different
available model with different scales and numerical costs. Implicit schemes
exist, but their suitability and speed here are not established by this audit.

### Gkeyll

Source inspected at `7c4b7b185d04d164afbe9aa8548331788646a538`; MIT license.
Its five-/ten-moment species models retain electron inertia, electromagnetic
coupling and (for ten moments) a pressure tensor. These are valuable independent
equations/reference implementations for the components under development.

The inspected standard [Euler primitive/flux functions](https://github.com/gkeyllorg/gkeyll/blob/7c4b7b185d04d164afbe9aa8548331788646a538/moments/zero/gkyl_wv_euler_priv.h)
and [ten-moment pressure functions](https://github.com/gkeyllorg/gkeyll/blob/7c4b7b185d04d164afbe9aa8548331788646a538/moments/zero/gkyl_wv_ten_moment_priv.h)
divide momentum products by density. We have not identified a ready-to-use
exact-vacuum/free-surface treatment for our problem in this path. That is a
specific integration gap, not a claim that every Gkeyll application needs a floor.

The inspected [moments application](https://github.com/gkeyllorg/gkeyll/blob/7c4b7b185d04d164afbe9aa8548331788646a538/moments/apps/moment.c)
builds its geometry on the host, and its Lua moment-equation constructors request
`use_gpu=false`. GPU support in the wider Gkeyll framework therefore does not
establish GPU acceleration of this particular fluid application. Its kinetic
branches are separate choices. No Gkeyll executable or runtime was benchmarked
in this bounded comparison; this is explicitly a source/capability audit.

### Numerical infrastructure and visualization

PETSc remains a candidate for accelerating our coupled sparse systems, subject
to profiling and a suitable block preconditioner. Installing it does not resolve
the missing physics. ParaView is suitable for field-line/volume inspection;
neither tool requires replacing our experiment definitions or physics API.

## Verification protocol

1. Common cold parallel Hall reference, n=1e18 m^-3, physical proton mass,
   B0=.001 T, perturbation=1e-7 T, k di=1 and both helicities. Periodic z,
   zero temperature/resistivity/hyper-resistivity, no filter, regular particle
   loading. Compare 32/64/128 cells at common final phase omega t=2.
2. Refined-step and particle-loading controls separate time/particle effects
   from spatial error. Compare the existing four-mode linear component against
   the same continuum solution; do not rank general runtime against that
   much smaller mathematical problem.
3. A full-Maxwell vacuum wave has no particle species or density floor. It
   travels one wavelength through a periodic domain. Compare stagger-aware
   continuum fields, refinement and endpoint electromagnetic energy.
4. Archive complete inputs, field samples, particle momenta/weights, source
   identities, energy diagnostics and timings. Recompute errors and energy
   from the stored physical arrays before making a recommendation.

Neither fixture tests a moving plasma-vacuum edge, cylindrical axes, collisions,
coil-driven formation or 3D stability. CPU/GPU speedup is not measured.

## Preserving expansion without committing to a replacement

The useful boundary is the experiment specification and result record. A future
WarpX adapter can translate the same prescribed sources while running in its
own process/environment. In our convention A_theta=psi/r, with its regular
axis limit, so the existing coil calculation supplies the vector potential
needed for an external-field adapter. Current sources are prescribed waveforms;
this does not introduce a hardware or circuit design.

Capabilities must be explicit: RZ versus 3D, quasistatic versus full Maxwell,
fluid versus kinetic particles, finite electron inertia, exact vacuum versus
regularized Ohm law, and each transport closure. Unsupported combinations must
reject rather than silently change density or physical assumptions. Output
conversion must retain native staggering, physical units, boundary conditions,
source work, model limitations and the engine/version identifier. Reconstructed
flux/topology must pass a divergence/consistency check before comparison.

There is no reason to discard the existing code to make that extension possible.
There is also no basis to assert that all remaining formation physics is a
small finishing task. The next engine decision should follow a common physical
plasma-edge/coil case, not familiarity or lines of code already written.

## Measurements and recommendation

**Retain the existing project and custom axisymmetric research path. Keep WarpX
as an optional independent/kinetic engine; do not replace the core now.**
This is a fit decision, not a claim that our incomplete formation solver is more
accurate in general. The clearest potential saving is future kinetic/3D particle
work. The exact-vacuum/edge and physical-regime gaps prevent calling either
candidate a verified shortcut to our present formation target.

All **17 final GPU cases complete and pass archive/physical replay**. Source
identities, full inputs, actual accepted time/step count, fields, particle
weights/momenta and native energy histories are retained in
`outputs/engine_comparison`. The figure was visually inspected.

Relative magnetic-field errors, expressed as percentages:

| Fixture | 32 cells | 64 cells | 128 cells |
| --- | ---: | ---: | ---: |
| Hybrid ion-cyclotron branch, fixed particle time step | 0.254055% | 0.151745% | 0.145373% |
| Hybrid whistler branch, fixed particle time step | 2.445204% | 0.647522% | 0.217297% |
| Empty full-Maxwell domain, fixed CFL | 0.756533% | 0.189219% | 0.047310% |

The Maxwell spatial orders are 1.99935 and 1.99984. At 128 cells the maximum
energy variation over the entire period is 1.134e-8 relative; the smaller
2.49e-15 endpoint change must not be mistaken for the maximum drift. Halving
the time step raises its continuum field error to 0.059138% while reducing
maximum energy variation to 7.085e-10. This is consistent with Yee spatial and
temporal dispersion partially cancelling: its dispersion relation is
sin(omega_num dt/2)=(c dt/dz) sin(k dz/2), so reducing the Courant number can
increase the net phase error at fixed mesh.

The hybrid results are **not a clean spatial-convergence result at the chosen
particle step**. At 128 cells, baseline/half/quarter-step successive field
differences have ratios 2.0009 (ion branch) and 1.9894 (whistler). This indicates
first-order temporal behavior in this tested configuration; its cause is not
established as an upstream defect. Correct continuum E initialization does not
remove it. Quarter-step continuum errors are 0.037975% and 0.162279%, respectively.
Perturbation energy increases by approximately 0.223%/0.112%/0.0559% as the
particle step is halved twice, independent of background energy dilution.
Native intermediate histories and independent particle/field endpoint energies
agree within explicitly recorded text-output precision. This is a limitation
to resolve before quantitative adoption, rather than a claim of exact conservation.

Increasing particles per cell 16→64 changes the finest whistler field error
from 0.002172967 to 0.002172318 (fractional units). Reducing perturbation amplitude
tenfold or increasing the inactive density floor by 1e6 leaves the normalized
error essentially unchanged. The floor control is entirely in positive-density
plasma and supplies **no evidence about vacuum correctness**.

The final small GPU evolutions take about 0.04–0.21 s for the vacuum cases and
1.6–16.9 s for the hybrid cases, including diagnostics and host readback;
initialization is separate. These tiny 1D runs do not estimate an FRC runtime,
3D capacity or CPU/GPU speedup. The GPU installation is isolated at
`/home/m/Code/frc-tool-comparison`, with a private `.venv`; no system packages
were changed. It remains available for future comparisons.

The existing four-mode Hall component passes the same continuum target with
relative state errors 1.588e-6→3.975e-7 (ion branch) and
1.084e-5→2.721e-6 (whistler) under time halving, with energy change around 1e-15.
Sixteen existing Hall-wave tests pass. This component has neither PIC particle
error nor a spatial grid, so those smaller errors and runtimes are **not** an
apples-to-apples engine performance ranking.

### Preserved setup controls

- `setup_trials_01` preserves two initial Hall runs and the vacuum startup
  rejection caused by omitted zero field components. Explicitly supplying all
  components fixes the deck; no physical vacuum evolution had occurred in the
  rejected run.
- `initial_electric_control` preserves the original 15-case suite, source and
  audit with an incomplete Hall initial electric field. WarpX desynchronizes
  velocity before the first hybrid current deposition, so initial E matters.
  The final suite supplies Ex=(omega/k) By, Ey=-(omega/k) Bx and Ez=0 from the
  analytic solution. We do not label initialization contamination a solver bug.

### Next action

Resume the existing finite-mass/current verification checkpoint: complete its
pending result audit and repair the stationary-bias precision failure before
moving integration. Keep the exact-vacuum requirement and formation gate.
Use this installed WarpX reference when evaluating the eventual kinetic/3D
path. Before adopting it for formation, a separate common RZ coil/plasma-edge
case must establish boundary physics, source work, temporal/spatial convergence
and transport validity. That case and any production backend adapter remain
unimplemented; no migration or performance saving is claimed yet.

### Reproduction

The tested environment and source hashes are in `environment.json`. From the
remote comparison root, the successful build used:

```sh
cmake -S warpx-26.09 -B build-1d-cuda \
  -DWarpX_DIMS=1 -DWarpX_COMPUTE=CUDA \
  -DCMAKE_CUDA_COMPILER=/opt/cuda/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=89 \
  -DWarpX_MPI=OFF -DWarpX_EB=OFF -DWarpX_QED=OFF \
  -DWarpX_OPENPMD=OFF -DWarpX_FFT=OFF -DWarpX_PYTHON=ON \
  -DWarpX_PYTHON_IPO=OFF -DCMAKE_BUILD_TYPE=Release \
  -DPython_EXECUTABLE=/home/m/Code/frc-tool-comparison/.venv/bin/python
cmake --build build-1d-cuda --target pip_install -j 4
```

The private environment needs pip for CMake's installation target and
`cupy-cuda13x` for Python GPU field access. Both are installed. In `bench`, set
`CUDA_VISIBLE_DEVICES` to the recorded RTX 4090 UUID, then run
`sh examples/run_warpx_comparison.sh` with a new `FRC_COMPARISON_RUN_ROOT` if
repeating the dataset. The launcher preserves existing records; archive replay
checks exact current driver/reference hashes. Locally run
`python -m examples.backend_custom_comparison` and
`python -m examples.backend_comparison_report` using the project's environment.
