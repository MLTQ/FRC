# 3D hybrid model of the formed FRC (frc-ibl.11)

## Why
The 2D formation solver cannot show the effects that set FRC lifetime or the
user's experiments around it:
- n=2 rotation and tilt;
- quadrupole stabilization;
- rotating-magnetic-field (RMF) current drive;
- particle confinement.

FRC ions in the 12" device are kinetic: s̄ ≈ 2–3, only a few gyroradii across
the separatrix. The standard self-consistent model for this regime is
**hybrid**: kinetic ions with fluid electrons. The first implementation starts
from a formed 2D FRC. The user's current target is full 3D evolution of every
configuration from the initial fill at 0 through 40 µs, including formation.

## Current measured state and next action
The original Option A continuation covers physical time 10–16.002 µs, on a
10×10×20 mm grid with 645,040 initial macro-ions. Its terminal n=0 flux is
1.979 mWb, radius 11.48 cm and length 78 cm; the nearest 2D state gives
1.846 mWb, 11.18 cm and 97.8 cm. Instantaneous midplane n=1/n=2 amplitudes
end at 1.31%/2.13%, without a clean exponential growth interval. Loaded-ion
loss is 31.98%, including initial gas outside the separatrix. Number plus
absorbed ions closes to 4.4e-16 relative.

This historical prototype predates the matched boundary-current support
described below. It is not a stability or confinement result. Ion kinetic + magnetic +
absorbed-ion energy rises by 1,075 J; the isothermal-electron and boundary
energy exchanges are not closed. M3 equilibrium and M4 stability remain open.
All 61 viewer frames reproduce saved times, contain finite current textures
and retain genuine 3D field/current samples. Ion dots have no persistent IDs;
the display-count slider selects a subset without changing the simulation.

All six corrected v4 shots completed actual 0–40 µs, with 40,000 ion steps
and 401 frames each. All 2,406 viewer states match raw snapshot times and
sample counts, with finite current textures, ion dots and field lines.
Evidence and exact source archives are in `outputs/hybrid3d_full_shots_v4/`.
The initial fill is rebuilt
at actual t=0; no evolved 0.1 µs snapshot is relabeled. The current driver
uses conservative electron internal energy, adjoint pressure compression,
Spitzer heating and counted wall/empty-node energy loss. No temperature,
energy or mode amplitude is rescaled. No azimuthal filter or mirrored sector
is used in evolution. The display's gold n=0 flux outline is diagnostic only.

**Measured correction.** An actual 4 µs startup with this closure has a
0.104 J load-energy accounting discrepancy, versus 419.8 J for the failed
entropy model. This includes measured cutoff removal, outgoing electron
enthalpy, energy left on emptied nodes and absorbed ions. It does not prove
those loss models are physically converged. Explicit pressure-centering
error and grid-equation residuals are recorded separately. A guard flags a
full shot above max(5 J, 1% of image work); it never repairs the state.

| Completed shot | Initial loaded ions remaining at 40 µs | Numerical heat removed | Load accounting residual |
| --- | ---: | ---: | ---: |
| Option A | 42.07% | 726.65 J | +0.348 J |
| End plugs | 50.51% | 731.50 J | +0.458 J |
| Re-boost | 41.06% | 779.50 J | +0.378 J |
| Inner tube | 0.766% | 526.91 J | +0.031 J |
| Flux core | 0.718% | 532.74 J | +0.032 J |
| Core punch | 0.713% | 528.03 J | −0.004 J |

These counts cover all initial gas, including material outside closed flux;
they are not a confined-cohort measurement. Particle number plus absorbed
number closes within 2.5e−14 relative. Numerical removals combine rejected
low-density Joule heat and energy on emptied PIC nodes: 43–60% of delivered
load work, a major unresolved model error. They are not physical radiation
or measured wall heat. The inner wall absorbs ions; the older 2D fluid wall
reflects them, so the two models cannot be compared as identical experiments.
No physical configuration ranking follows from this table yet.

The table preserves raw v4 histories and separately accounts for the frozen
reference flux of banks connected after formation. End plugs have a raw
image-work residual −9.376 J plus 9.834 J reference exchange. The posthoc
0.25 µs reference quadrature leaves source waveform lag in the residual.
An independent actual 0–4 µs end-plug replay with direct ion-step source work
gives +0.1065 J error and bitwise-identical B, n and pe to the original v4
states. This is an accounting change only, not an energy correction to the
simulation. See `source_work_verification.json`, per-run
`energy_reference_audit.json`, and `summary.json` in the output campaign.

**Measured sensitivity, not convergence.** At 0.75 µs, halving spacing
(20×20×40 to 10×10×20 mm) and ion timestep (1 to 0.5 ns), at eight ions per
cell and 5% cutoff, changes response magnetic energy 344→466 J, ion kinetic
energy 197→117 J and numerical removals 20→106 J. On the coarse grid, reducing
the cutoff from 5% to 1% gives 329 J, 195 J and 10 J respectively. Formation
energy partition is unresolved even though the energy ledger closes.
`startup_sensitivity.json` records the actual runs. Next: replace or justify
the numerical losses and measure formation grid/time/cutoff/particle-count
convergence, then evaluate M3/M4 against independent equilibria/growth rates.

**Requested higher-resolution campaign.** `examples/hybrid3d_campaign.py`
now runs all six actual 0–40 µs shots at 10×10×10 mm, eight ions/cell and
0.5 ns ion steps. The baseline grid is 39×39×232 including padding, with
161,260 gas cells and 1,290,080 initial macro-ions, versus 9,900 gas cells
and 79,200 ions in the coarse campaign. Transverse spacing is halved and
axial spacing quartered; the ion timestep is halved, with the field timestep
still determined by its existing stability limits. This changes multiple
resolution parameters together and is not an isolated spatial-order proof.
Physical inputs, 5% density cutoff and output interval 0.1 µs are retained;
no symmetry constraint or artificial disturbance is added.

The first CPU queue, `outputs/hybrid3d_shots_r10/`, failed at 1.141 µs:
outgoing electron-energy transport made an empty near-wall node negative.
The field whistler/acoustic timestep did not bound multidimensional advection,
wall enthalpy and pressure expansion on partial gas volumes. Its other workers
were stopped, and its incomplete evidence is retained. Electron transport now
uses its own conservative positivity substeps; see the derivation below and
`plasma_bench/hybrid3d_energy_electrons.md`. This numerical change means the
new comparison is not a pure grid-convergence experiment against v4.

The corrected six-shot queue is `outputs/hybrid3d_shots_r10_cfl_metal/`
(`gpu=2`), with a matched CPU64 baseline in
`outputs/hybrid3d_shots_r10_cfl_cpu_control/`. Both started from actual fill,
target 40 µs and import immutable source copies with SHA256 manifests.
The Metal queue permits two workers, two Numba threads each, but releases
other cases only after baseline reaches 1.3 µs, beyond the old failure.
A baseline failure defers queued cases. The CPU control uses one worker.
No full corrected high-resolution shot has yet completed. A load accounting
guard invalidates a failed shot without correcting its state. Queued shots
defer below 10 GB free disk. Completed shots automatically join the existing
viewer; coarse runs remain. `campaign.json` records actual saved times and
failures; `DONE` only means that the queue has reached terminal states.

**Why the current can still look smooth.** The older fluid viewer rotates
an axisymmetric solution. The hybrid current volume instead uses the real
Cartesian response field. At 40 µs in Option A, a fixed interior cylindrical
sample has nonaxisymmetric current L2 fraction 0.363 (13.2% of squared norm),
with m=1/m=0 RMS 0.168 and m=2/m=0 0.124. This includes particle/grid noise;
the m=4 component can include Cartesian mesh bias and is not instability
evidence. Shared color scaling and the gold n=0 outline conceal some of
this variation. A current-volume screenshot alone does not establish magnetic
topology. Published 3D FRC calculations can develop nonlinear tilt and
rotational modes ([Belova et al., PPPL-4075](https://bp-pub.pppl.gov/pub_report/2005/PPPL-4075.pdf));
the user's unidentified reference image cannot be assigned a particular mode,
formation history or physical model. We do not seed a distortion to match it.

**Independent thermal comparisons, production default.** A moving-density,
uniform-pressure contact translates one period with pressure error <4e−16
and no net energy change at16/32/64 cells. Force-free constant-η decay errors
are0.64%→0.16%→0.040%, and Spitzer errors0.49%→0.12%→0.029%, against an
independent continuum ODE. Total-energy error is below5.3e−9 (constant η)
and3e−12 (Spitzer), relative to initial magnetic energy; ions remain at rest.
Evidence: `outputs/hybrid3d_energy_check_production.json`. These periodic
runs do not verify formation or an absorbing boundary.

**Prior failures, preserved.** The v1–v3 campaigns are invalidated. Per-ampere
coil interpolation lost source amplitude; unresolved coil curl entered Ohm's
law; conductor currents on excluded E edges injected Hall energy; entropy
transport then manufactured heat. The last defect accounts for419.7 J of
419.8 J full-shot error at4µs (`outputs/hybrid3d_pressure_work_diagnostic`).
A positive entropy transfer did not cure the incompatible pressure work.
Those runs and the first positivity-failing grid-energy prototype must not
be interpreted as physical configuration comparisons.

**Vacuum, geometry and drive.** An absorbing inner tube contains exact
vacuum and retains the meridional circulation lost by a single-valued
potential solve. A current-loop field comparison converges approximately
at second order. Known coil bases are discretely current-free in gas,
retaining boundary-normal flux and compatible A; the correction changes
source B by0.894%→0.314% at20→10mm. Only the known source is constrained.
Plasma current and E use matching edge support and gas-volume node weights.
Projected-coil mutual linkage errors against independent analytic filaments
are0.527%→0.151% (solid) and0.498%→0.140% (annular) at20→10mm. The switched
banks use the same flux-boundary/image-load model as the 2D formation path,
now driven by full Cartesian plasma current; this is not an experimental
validation of the multi-coil 3D boundary model.

**Remaining limits.** No full-shot grid/cutoff convergence or 3D instability
benchmark is complete. Measured numerical heat removal at the cutoff and
empty nodes dominates the uncertainty and is counted separately from ion absorption.
Upwind electron-energy transport and explicit time splitting are first
order and diffusive. CIC, node/edge averaging, RK3 and staircase boundaries
add error/diffusion even with η_H=0. Spitzer uses lnΛ=10 and max(Te,1eV);
pressure has no temperature floor. Ionization, electron conduction,
temperature equilibration and sheath physics are absent. Closed n=0 flux
and a bounded energy budget do not establish 3D stability or confinement.
Legacy isothermal/barotropic choices remain available; the entropy closure
is retained only for reproduction (`closure='entropy'`).

## Model
**Ions.** Full-orbit macro-particles (H⁺ or D⁺), advanced by the Boris push
with trilinear (CIC) weighting. Density n and flux n u_i are deposited on the
grid nodes.

**Electrons.** A massless, quasi-neutral fluid (n_e = n_i) obeying the
generalized Ohm's law, with the Hall term and electron pressure:

    E = −u_i×B + J×B/(n e) − ∇p_e/(n e) + η J − η_H ∇²J,   J = ∇×B/μ0

- Thermal default: conservative internal energy U_e=p_e/(γ−1), with
  ∂tU_e+div[(U_e+p_e)u_e]=u_e·grad(p_e)+ηJ² and u_e=u_i−J/(ne).
  Legacy isothermal/barotropic choices are explicit alternatives.
- η is physical resistivity (uniform or temperature-dependent Spitzer).
- η_H is an optional **declared** hyper-resistivity that damps grid modes.
  It is zero in the current campaign; the thermal model accounts for ηJ².

**Fields.** Darwin, with no displacement current.
- Faraday's law, ∂B/∂t = −∇×E, on a Yee grid: B on faces, E and J on edges.
- ∇·B is preserved to round-off.

**Low density.** Ohm's law uses max(n, n_floor). The floor is declared and its
heat not retained below the cutoff is counted. A vacuum region is not a
plasma; see the geometry milestone.

## Time integration
Per ion step dt:
1. E^n is computed from Ohm's law with B^n, n^n and u^n. Here
   u^n = 1.5 u^{n+½−1} − 0.5 u^{n−½−1} is extrapolated from the two most
   recent half-step fluxes.
2. Boris push: v^{n−½} → v^{n+½} with E^n and B^n at x^n, then
   x^{n+1} = x^n + dt v^{n+½}.
3. Deposit n^{n+1} at x^{n+1}, and the flux at x^{n+½} with v^{n+½}.
4. Advance B from t^n to t^{n+1} in subcycles of SSP-RK3. Ohm's law is
   re-evaluated every stage with n interpolated in time and u = u^{n+½}.
   RK3's stability on the imaginary axis (|λ dt| ≤ √3) covers the dispersive
   whistler.

**Time-step limits:**
- ions: dt ≤ 0.1/Ω_i, and a particle travels < 0.5 cell per step;
- fields: dt_B ≤ √3 / ω_max, with ω_max = v_A k_max (1 + k_max d_i) at the
  grid scale.

## Geometry (milestones)
- **M1:** periodic box, for verification.
- **M2:** the 12" device, on a Cartesian grid, which avoids an axis
  singularity.
  - The glass tube (r_w = 15.24 cm) and its ends absorb ions.
  - The single-turn coil is a perfect conductor at r_c = 16.5 cm, with
    conducting planes beyond the tube ends, as in the 2D model.
  - Exact vacuum fills the gap between tube and coil (no ions; J = 0). How it
    is enforced is decided in M2 and declared. Any approximation is measured
    against the 2D solver's exact vacuum.
- **M3:** start from the formed FRC of the 2D solver (flux, density, pressure,
  velocity) and from an analytic rigid-rotor equilibrium. Check that the n=0
  equilibrium holds, and compare n=0 evolution with the 2D run.
- **M4:** n=2 rotation and tilt; a quadrupole at 7% of B_e.
- **M5:** RMF antennas, even and odd parity.
- **M6:** hybrid runs and ion trails in the 3D viewer.

## Verification (independent, no unit tests)
- **Particle orbits:** the Boris gyration phase error against the exact orbit,
  and exact energy in static fields.
- **Parallel waves** in a uniform cold plasma, against the exact Hall-MHD
  (cold-ion hybrid) dispersion:

      ω² = k² v_A² (1 + ω/Ω_i)

  This gives the whistler branch (ω > 0) and the ion-cyclotron branch
  (ω → −Ω_i). Measure the frequency for k d_i = 0.25–2 and its convergence
  with the grid.
- **Conservation:**
  - ∇·B to round-off;
  - total energy (ion kinetic + magnetic, cold electrons, η = 0) against step
    size.
- **Later:**
  - the rigid-rotor equilibrium (M3);
  - published n=2 and tilt growth rates (M4);
  - RMF penetration and current against theory (M5).

## Performance
Particle push/deposit and device field kernels are JIT-compiled (numba).
The completed coarse six-run campaign used one thread per run; the new
10 mm Metal queue allows two shots, two threads each, plus one CPU control.
Its startup gate initially permits only the baseline. This machine is
an Apple M1 Pro with eight performance/two efficiency CPU cores, a 16-core
GPU and 32 GB unified RAM. One observed startup snapshot used about 110% CPU
and 0.8 GB resident memory per shot: configured threads are not fully utilized
throughout every step.

**Measured CPU scaling.** Actual identical 10 mm initial-fill shots were
replayed to 0.02 µs with 1/2/4 threads while the two full shots continued.
The timed interval 0.005–0.020 µs excludes initialization and ten warmup steps;
each advances 30 ion steps and 120 field substeps, with cProfile enabled.

| Threads per shot | Measured seconds | Throughput versus one thread |
| --- | ---: | ---: |
| 1 | 23.185 | 1.000× |
| 2 | 19.178 | 1.209× |
| 4 | 18.437 | 1.258× |

Four threads improve throughput only 4.0% over two in this short startup
measurement. Response magnetic, kinetic and electron energies agree to
roundoff; no physical parameter or precision changed. This single pass with
background work is not a formed-plasma or total six-shot throughput benchmark.
Evidence: `outputs/hybrid3d_parallel_check/comparison.json` and `.prof` files.
In the one-thread profile, 40 full steps take 30.8 s: field `rate` calls
8.7 s, three particle-moment deposits per step 5.2 s, electron `advance`
4.0 s and vacuum `project` 3.1 s. These cumulative times overlap inside
their parents and must not be summed with the parent step time. The sparse
LU solve itself takes 1.45 s. Array allocation, reductions, energy transport
and boundary work limit scaling beyond the existing parallel kernels.

**Live profile during the fine run.** Brief native stack samples of the
two active baselines around 0.8 µs leave both shots running. In the Metal
main thread, 168/232 samples (72.4%) are in NumPy arithmetic, indexing or
copying; 41/232 (17.7%) are compiled CPU work/worker synchronization,
9/232 (3.9%) MLX API/synchronization and 10/232 (4.3%) the sparse LU call.
The CPU control has 181/241 (75.1%) in NumPy. These single three-second
samples are indicative, not precise Python function attribution or full-shot
timing. Worker waits include useful work on other threads. Vacuum projection
also includes array work outside its sparse LU call. The next performance
candidate is reducing temporary arrays and fusing field/energy calculations;
porting only the sparse solve cannot remove most of the observed cost.
Evidence: `outputs/hybrid3d_live_profile/` (raw samples and classified summary).

**Float64 array fusion (frc-ibl.11.12).** Optional `arrays=1 thermal=1`
(`--arrays` on a new campaign) fuses the existing compression, upwind flux,
pressure-power, edge-resistivity average, Joule transfer and energy RHS loops.
Scratch storage removes intermediate arrays; pressure-work and heat values
retained by RK stages remain separately owned. There is no fast math,
parallel reduction, equation change, symmetry constraint or physical rescale.
Spitzer evaluation, all budget/positivity checks, geometry, exact vacuum,
particle integration and circuits retain their previous implementation.

| Matched replay, two CPU threads | Original | Fused float64 arrays |
| --- | ---: | ---: |
| CPU, 20×20×40 mm, 0–0.75 µs, 1,490 warmed steps, cProfile | 66.16 s | 57.27 s |
| CPU, 10 mm, 1,290,080 ions, 50 warmed startup steps | 45.86 s; repeat 38.67 s | 31.47 s |
| Metal+CIC, 20×20×40 mm, 0–0.75 µs, 1,490 warmed steps | 96.24 s | 78.08 s |
| Metal+CIC, 10 mm, 1,290,080 ions, 50 warmed startup steps | 35.25 s | 28.39 s |

The fine warmed interval is 0.005–0.030 µs. Initialization and the first ten
steps/JIT are excluded; all remaining work and synchronization are included.
CPU fine reference measurements bracket the candidate: it uses 18.6–31.4%
less time; Metal fine uses 19.5% less (1.241× throughput). This is a busy
shared host with varying load, not a full-shot or guaranteed speedup.
An initial fused-loop version had no benefit; conditional periodic neighbors
and explicit component unrolling removed expensive inner-loop indexing.
The measured compiler executes the unrolled CPU transport/Joule loops
serially despite their parallel hint (performance warnings are retained in
the logs); the other four kernels parallelize. This is included in the
measurements, not a claim that every array pass uses both threads.

Complete CPU64 B, all ion positions/velocities/weights, n, U and both flux
arrays are exactly equal to the pre-change implementation at both endpoints.
The 0.75 µs load residual remains −0.104644578 J; roundoff correction is zero.
Existing independent Gaussian/contact/constant-eta/Spitzer comparisons pass
or converge with unchanged errors. Final Metal pair differences are ≤1.17e−7
in total B relative L2, 1.31e−7 in n and 9.18e−8 in U at 0.75 µs, compatible
with the existing atomic deposition precision; the long-time sensitivity
below remains unresolved. Default stays off and current source-frozen queues
are not restarted. Evidence: `outputs/hybrid3d_array_benchmark/`, including
the pre-change source, complete float64 states, comparisons and final overlay.

**Working Metal prototype (frc-ibl.11.11).** MLX/mlx-metal 0.32.4 are installed
and pinned in the optional `metal` extra. `gpu=1 thermal=1` selects fused GPU
transport fluxes and wall rates; `gpu=2` also selects particle deposition.
GPU arithmetic in these blocks is float32. Particle positions/velocities,
stored electron energy, conservative divergence, compatible pressure work,
Faraday curl, exact vacuum solve and circuits remain CPU64. CPU remains the
default; the corrected six-shot queue explicitly selects GPU mode 2. The new modules
are `hybrid3d_metal.py` and `hybrid3d_metal_energy.py`.

The deposition kernel scales moment values, combines matching destinations
within eight-lane groups, and falls back to general atomics for other keys.
This changes arithmetic scheduling, not particle ordering or geometry.
An initial uncombined atomic implementation was slower than CPU; the final
combining implementation is measured below. Speed depends on particle
locality, and results are not an estimate of late-shot or full-shot speed.

| 10 mm device startup, two CPU threads | CPU64 | Metal transport + CIC |
| --- | ---: | ---: |
| Complete 30 steps, 0.005–0.020 µs | 16.57 s | 14.42 s |
| Actual 1,290,080-ion deposition, including conversions/sync | 21.59 ms | 6.79 ms |
| Deposition with inputs resident, kernel/sync only | 21.59 ms | 1.80 ms |

Complete-step throughput improves 1.149× in this startup sample (13% less
time). Deposition improves 3.18× including movement and 12.00× for resident-input kernel work
alone. Those kernel ratios are not full-solver ratios. Initialization and the
first ten JIT warmup steps are excluded; all remaining CPU fields, boundary,
circuit work and GPU conversions/synchronization are included in step timing.
Transport-only takes 16.69 s, slightly slower than CPU. These are matched
measurements after the transport-CFL correction. Background load changed
between earlier prototype versions, so their absolute times are not a kernel
revision speed comparison. Evidence: `outputs/hybrid3d_gpu_cfl/`.

**Physical comparisons.** Both modes completed moving-contact and constant/
Spitzer force-free heating runs at 16/32/64 cells. Transport-only contact
pressure error is ≤2.1e−15; adding float32 CIC gives ≤3.9e−8. Magnetic decay
errors decrease 0.641%→0.161%→0.0402% (constant η) and
0.490%→0.120%→0.0290% (Spitzer), against independent continuum solutions.
Energy errors stay below 5.3e−9/2.6e−12 of initial magnetic energy.
Independent ordered/shuffled uniform-current deposits reproduce exact node
density/current; a partial final SIMD group closes number to 1.7e−16.

**Independent transport timestep.** For gas dual-volume fraction f, upwind
outgoing rate A, wall rate b and central volume-velocity divergence D(fu),
the own-node energy-loss rate is
λ=[A+γb+(γ−1)max(D(fu),0)]/f. The shared CPU/GPU update uses
ceil(dt·max(λ)/0.8) equal electron-energy substeps. This leaves a nonnegative
own-node coefficient, incoming flux and heat. Velocity and stage heating are
fixed during these split substeps; pressure evolves and its centering error
is measured. There is no energy/temperature rescale or material-negative
clamp. This is first-order upwind transport, with numerical diffusion; it
does not fix the existing cutoff/empty-node heat sinks.

A periodic Gaussian translation at outer transport CFL 1.45–1.49 requires
two energy substeps. CPU and Metal at 32/64/128 cells remain positive, with
no roundoff energy correction and total-energy relative error ≤2.9e−16.
Relative L1 translation error decreases 0.825→0.591→0.387: diffusion is large
on this narrow Gaussian but decreases with refinement. This prescribed-flow
operator fixture is not a force-balanced formation experiment. Evidence:
`outputs/hybrid3d_gpu_cfl/energy_advection.json`.

Actual coupled-bank compression replays reached 0.75 µs on 20×20×40 mm,
0.5 ns ion steps. Endpoint load-budget errors are −0.104644578 J (CPU),
−0.104644578 J (transport) and −0.104638516 J (transport+CIC). The latter
differs from CPU by at most 1.024e−5 in saved B-component relative L2,
4.64e−7 in n and 2.95e−6 in p_e; scalar energy differences are below 8 µJ.
Saved arrays are float32 and comparisons use float64 norms. These coarser
runs take 45.96/45.00/45.87 s respectively: there is essentially no complete
step benefit at this smaller grid. All three require at most two energy
substeps, with maximum actual loss fraction 0.78221 and zero energy roundoff
correction. CPU pressure-centering residual is 0.00570 J and the grid-energy
identity residual is −1.31e−13 J. Evidence:
`outputs/hybrid3d_gpu_cfl_compression/`, including raw endpoints,
comparison summary and exact source hashes/archive.

The corrected 10 mm full-shot controls have also reached their first saved
0.1 µs (200 steps). CPU/GPU load residuals are −0.003116730/−0.003116698 J
(difference 3.22e−8 J), with no electron roundoff correction. Saved relative
L2 differences are ≤8.35e−6 for B components, 1.90e−6 for n and 1.96e−5
for p_e. Loaded-ion plus absorbed-ion count balances to 2.11e−16 in each
history; this count is separate from the atomic deposition error measured
above. First-interval wall times are 125.31/119.45 s, about 1.049× GPU
throughput, including cold step/JIT work with concurrent shots. This is not
a full-shot timing. Evidence:
`outputs/hybrid3d_shots_r10_cfl_metal/matched_cpu_0p1us.json`.

**Unresolved mixed-precision sensitivity in the fine run.** The short
comparisons above do not establish long-time local-field agreement. Saved
CPU/Metal response-face relative L2 differences at 0.1/0.5/0.8/1.0 µs grow
4.87e−6→0.00197→0.01069→0.09501. At 1 µs, response-node difference is
0.04391, density 0.000886 and pressure 0.03052. Using the mean CPU/GPU
density for unweighted node regions, field differences are 0.002244 above
the declared cutoff, 0.13712 at nonempty nodes below it and 0.04369 at empty
nodes. These regional norms are diagnostics, not physical energy integrals.
At 0.8 µs the scalar energy differences are still below 0.003 J and load
accounting residuals differ by only 0.000142 J. Thus energy-ledger agreement
does not establish agreement of detailed 3D structure. These measurements
do not identify whether precision, amplified noise, timestep sensitivity or
an unresolved numerical mode dominates. The full CPU control continues;
the GPU configuration campaign remains exploratory. Investigate the
low-density amplification and physical heat-removal treatment before using
fine current topology or configuration rankings as predictions. Evidence:
`outputs/hybrid3d_shots_r10_cfl_metal/mixed_precision_drift.json` and
`matched_cpu_0p8us.json`. Claude owns the physical-loss follow-up on
frc-ibl.11.9; Codex retains performance/run coordination.

The same controls now have a matched 3.9 µs comparison. Response-face L2
difference is 0.2301, above-cutoff nodal field difference 0.1603, density
0.1589 and pressure 0.2257. Global magnetic energies differ by 1.54 J
(0.14%), while load-accounting residuals remain +0.0249/+0.0099 J and
numerical heat removal is 547.7/546.7 J against about 1,498 J of source work.
The sensitivity has therefore reached the denser region as well; it cannot
be dismissed as an empty-node rendering artifact. These single-seed paths
still do not identify which trajectory is accurate or isolate noise from
precision/timestep effects. Evidence:
`outputs/hybrid3d_shots_r10_cfl_metal/mixed_precision_drift_3p9us.json`.
The 2026-10-22 00:02 EDT status check extends the concern to bulk n=0 geometry:
at matched 10.8 µs, Metal/CPU trapped flux is 1.3765/0.9009 mWb (52.8% higher
relative to CPU), and length is 0.72/0.54 m. Numerical heat removal is nearly
identical at 746.10/746.39 J. Neither path is established as accurate. Latest
saved progress is baseline Metal 11.6, CPU 10.8 and end plugs 7.8 µs of 40; four
Metal configurations remain queued. Baseline's latest 20 intervals give 0.937
simulated µs/hour, extrapolating about 30 hours remaining at unchanged load;
this is not a completion guarantee. Ledger residuals remain below 0.44 J, but
numerical removal is 761.7/1481.9 J (51.4% of source work) on the baseline.
Evidence, exact histories and matched rows: `outputs/hybrid3d_status_20261022/status.json`.
Partial saved-frame exports are available in the viewer's Run menu for both
controls and end plugs, with their actual time spans and 40 µs target. Metal
entries explicitly disclose their mixed-precision development status. These
exports are snapshots of progress, not an automatically refreshing stream.

**Backend isolation and coordinate prototype (2026-10-22).** Five actual fine
0–0.1 µs replays use identical 1.29M-ion initial fill, 10 mm cells, 0.5 ns,
coupled banks and unchanged loss physics. CPU repeat is strictly bitwise
identical in all native fields/particles/fluxes; its float32 B/n/pe also match
the existing frozen long CPU control exactly at 0.1 µs. At that endpoint,
response-face relative L2 differences from CPU64 are 7.16e-11 for GPU transport,
4.869e-6 for GPU deposition alone and 4.869e-6 for both. Early backend drift is
therefore dominated by deposition in this window, not by transport arithmetic.

The isolated `SplitCoordinateMetal` prototype computes particle cell/fraction
in CPU64, then sends integer cells and float32 local fractions to otherwise
unchanged float32 weights/velocities/atomics. It avoids subtracting global
positions/origins in float32. Independent number, first-position-moment,
linear-electrostatic-quadrupole work and uniform/shuffled/partial-group
comparisons were measured on actual native particles. Raw evolved-particle
density/flux L2 differences fall about 15/53x; atomic nonrepeatability remains
about 8e-8. Global manufactured-power error is not improved (~2e-9 of absolute
particle power), and no whole-step speed claim follows from these timings.

Two further actual coupled replays complete 0–0.1 µs with revised deposition
alone and combined with GPU transport. Response differences are 5.605e-9 and
5.613e-9 (~870x smaller than original combined); density is 1.08e-7 and pressure
9.84e-9. Source work remains 24.598 J, ledger residual about -0.003117 J and
numerical heat removal 0.021813 J. No energy/temperature/mode is fitted. The
prototype stays separate from normal production imports and existing queues.
A source-frozen combined prototype replay completes actual 1 µs (2,000 ion
steps) using the existing CPU control's matched fields. At 0.5 µs, original
and revised response-face relative L2 differences are 0.00196788 and
0.0000321292 (61.25x reduction); at 1 µs they are 0.0950132 and 0.0141901
(6.70x reduction). The earlier 870x factor does not persist. Revised 1 µs
dense/thin/empty nodal field differences are 0.0001089/0.0152155/0.0072004;
density and pressure differences are 0.000024398/0.00319827. CPU and original
GPU references are saved float32, revised states native float64. All three
native states are finite and source hashes intact. Source work is 1,001.529 J,
ledger residual -0.26694 J and cutoff-plus-empty-node removal 157.169 J
(CPU 157.298 J). No electron roundoff energy is added. Arithmetic agreement
improves while the physical-loss problem remains. A same-seed corrected repeat
is the next bounded action to separate atomic sensitivity from other errors
before any longer adoption. This does not establish physical loss/confinement
accuracy. Evidence:
`outputs/hybrid3d_precision_isolation_fine_v1/`,
`outputs/hybrid3d_cic_coordinates_endpoint_v2/`,
`outputs/hybrid3d_precision_split_fine_v2/`; the follow-up is
`outputs/hybrid3d_precision_split_to1us_v3/`.

Claude's fixed-mesh ion-count study uses actual counts 8, 27 and 64; earlier
nominal 16/32 requests both loaded 27. The loader now rejects non-cubic
counts and records the lattice/sample count without changing valid loads.
His 1 µs numerical removal decreases 68→53→45 J, and n=0 geometry agrees
closely, but the full 3D state is not established as converged. A read-only
comparison of his saved fields gives response-node differences 14.1→12.0%
above the cutoff and 52.2→51.4% below it for 8→27→64. These unweighted
pairwise norms use different thermal samples, not a paired ensemble proof.
They are in `outputs/hybrid3d_loading_contract/particle_count_field_sensitivity.json`.
See `examples/vacuum_mask_replay.md` for the physical-loss diagnosis; the
failed mask remains unadopted and the physical edge treatment remains open.

Claude's subsequent carrier ledger keeps fields unchanged: orphaned edges with
no ions in their four neighbouring cells dissipate 43.1 of 95.4 J of startup
Ohmic heat by 1 µs, and 111 of 113 J in the analytic compression. The isolated
electron-edge transport block lowers empty-node energy removal but raises
orphaned Ohmic dissipation: total startup removal 68.0→74.2 J, versus
182→144 J in compression. It is not adopted. The empty-node resistivity depends
on leaked energy before node cleanup rather than a physical electron closure.
His completed 5/2/1% floor/ledger comparison at matched actual 1 µs gives
cutoff-plus-empty-node removal 68.05/48.37/39.47 J (1% seed 1: 40.54 J).
Orphaned-edge Ohmic heat falls 43.1/21.0/12.2 J, while empty-node removal
rises 20.2/26.3/27.6 J and boundary electron enthalpy outflow 38.7/55.5/73.6 J.
Field substeps rise 3508/7510/14186. The midplane radius agrees between
2% and 1% at 0.112 m, versus 0.121 m at 5%; neither the other energy channels
nor the full 3D state are thereby established as converged. These coarse-grid
startup comparisons do not select the correct closure for the fine 40 µs runs.
The next bounded request is 1% plus the electron-edge transport block, against
the saved 1% controls. The passive ledger's sum of |J| times edge-dual area is
an extensive grid statistic, not a coil or net current; its stored values are
amperes (the orphaned values 8.36e6/7.63e6/6.61e6 A are MA, not kA).
The floor sensitivity remains a numerical-closure comparison, not an automatic
vacuum-interface solution. See
`plasma_bench/hybrid3d_carrier_ledger.md` and local task frc-ibl.11.9.

Full 0–40 µs GPU performance/precision and late-particle locality remain
unmeasured. The original numerical heat-removal uncertainty is still present.
A larger gain requires a resident particle/field/RK pipeline rather than
repeated host/device transfers, with compatible work and exact vacuum retained.
MLX GPU float64 support is unavailable; the precision split is explicit
([MLX data types](https://ml-explore.github.io/mlx/build/html/python/data_types.html)).

## Code

Oct 25 queue audit: both 10 mm baseline controls completed actual 40 µs, but
the main coordinator and end-plug/re-boost workers are absent. Last saved times
are 17.2 and 5.2 µs, respectively, on Oct 22; final native snapshots are finite
and no traceback establishes a cause. The three later configurations never
started. Stale running labels were corrected without launching/stopping jobs.
The snapshot's 4,000 sampled ions are for visualization, not a complete kinetic,
RNG or circuit checkpoint; do not claim an exact restart from these files.
Evidence and original manifest: `outputs/hybrid3d_status_20261025/`.

- `plasma_bench/hybrid3d.py`: grid, Yee operators, Ohm's law, Boris push,
  deposit, stepper.
- `examples/hybrid3d_verification.py`: the checks above.
