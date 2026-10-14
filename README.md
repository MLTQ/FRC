# Virtual plasma bench

The latest simulation couples the moving plasma to a magnetic-only exterior
and lets its enclosed flux evolve. Five complete pulse runs include spatial/time
refinement and evolving resistance. At the source peak, the plasma mean field
reaches 52–58% of the applied field; evolving resistance changes the heating by
+8.8% in this configuration. Read the
[measured results](docs/VACUUM_CYLINDER_RESULTS.md),
[response plot](outputs/vacuum_pulse_comparison/vacuum_coupled_response.png), and
[current handoff](docs/NEXT_SESSION.md). Finite coils and physical axial ends
are next; the present source and axial boundaries are idealized.

The preceding full runs establish [moving thermal relaxation](docs/ADVECTED_COULOMB_RESULTS.md)
and [magnetic decay with evolving resistance](docs/SPITZER_DECAY_RESULTS.md)
within their stated model limits. The sandbox is still a research prototype;
exact moving vacuum coupling, physical coil/end boundaries and credible FRC
formation remain unfinished. It is not yet validated for choosing a laboratory
formation experiment. Work uses actual simulations and mathematical baselines;
no unit tests are added or run under the current preference.

The [interactive 3D viewer](outputs/current_state_display_v3/index.html) still
shows the earlier nonuniform contact through 28.33 µs. It is a snapshot of that
run, not the newer driven-pulse calculation. Read its
[display limitations](docs/CURRENT_STATE_DISPLAY_V3.md).

The sections below preserve earlier component history; the current scientific
boundary is in [GOAL_STATUS.md](docs/GOAL_STATUS.md).

A scriptable, verified axisymmetric research sandbox. It currently implements
analytic vacuum electromagnetics and stationary resistive plasma induction with
a coupled nonconducting exterior. A developmental conservative MHD kernel now
simulates compression, axial motion, pressure/energy evolution and resistive
heating in explicit verification domains. The coil-driven formation workflow
now has a separate fixed-chamber coupling reference with moving gas, prescribed
coils and an elliptic magnetic vacuum. A new material-mesh reference evolves a
free plasma surface against a magnetic vacuum extending to infinity, with
scalar resistive induction and positive local Ohmic heating. Surface-focused
meshes now expose substantial spatial and time sensitivity and a conditional
thermal incompatibility of the initialized heated edge. Spatial coil coupling,
physical surface/transport models and appropriate end conditions remain
unfinished. Stability and particle confinement have not been established.
A separate electron-current component now includes finite carrier inertia and
its kinetic energy, with independent decay/screening references. Conditional
Coulomb and local Hall/inertia diagnostics are available. A bounded nonlinear
compressible Hall solver now evolves toroidal fields, rotation, density and
energy with compatible magnetic work and implicit Hall/resistive induction.
Separate electron/ion temperatures now evolve with moving ions in an ideal
reference and an experimental warm Hall coupling. Local enthalpy/Poynting
transport and fast warm waves refine, but finer tests expose a pressure-split
failure and substantial slow-wave error. Physical heating/transport and free-
edge coupling remain unfinished; this is not yet an accurate formation solver.

The [finite-mass current/thermal study](docs/CURRENT_THERMAL_VERIFICATION.md)
now audits 41 baseline trajectories and 33 precision controls. Physical electron
mass improves the original bounded pressure fixture; a stationary-background
roundoff rejection is fixed without suppressing weak waves. This component
still holds density and bulk flow fixed. The new
[moving-current study](docs/MOVING_CURRENT_VERIFICATION.md) verifies the exact
finite-mass momentum/energy reduction and 19 varying-density induction-rate
cases. The [conservative moving study](docs/FINITE_MASS_VERIFICATION.md) now
couples density, momentum, both species temperatures and full finite-mass
induction through an unsplit update. All 21 rate tests and 18 bounded trajectories
were audited, with 328 passing tests at that checkpoint. Global conservation is near roundoff, but
small ion-temperature changes and slow compressible modes remain inaccurate.
The measured reconstruction error, time controls and higher-order diagnostic
are retained; the physical plasma-vacuum edge and formation gates stay open.

The [smooth transport study](docs/SMOOTH_FINITE_MASS_VERIFICATION.md) now audits
31 rate cases and 27 trajectories, with 336 passing tests at that checkpoint. Shared transport
reduces fast-wave ion-pressure error from about 60% to 5.4% at the same grid.
Independent full-transit contacts include initial projection, continuous
temperature, time-step and flow-direction checks. Slow magnetic response still
fails because the velocity/current representations lose a delicate bulk/Hall
cancellation. These results do not establish accurate formation.

The [compatible field study](docs/COMPATIBLE_FIELD_VERIFICATION.md) now audits
25 improved linear cases and 25 preserved precision-failure comparisons, with
343 passing tests at that checkpoint. Higher-order magnetic/current fields and exact separate
bias storage reduce the finest slow-mode magnetic-rate error to about 0.00024%.
All three wave branches refine, with independent current/electric/energy checks.
This component is not yet connected to the moving plasma solver.

The [changing-density study](docs/VARIABLE_CURRENT_VERIFICATION.md) now extends
it with physical variable-density current inertia and local energy derivatives,
using the actual continuity rate. Its checkpoint passed 351 tests. Fourteen manufactured
component cases show strong magnetic/current-energy convergence. The corrected
species pressure-work pairing conserves its stated lumped energy, but the finest
entropy-rate error remains about 6.35%, concentrated near the axis. A consistent
higher-order thermal representation was the next target.

The [physical thermal study](docs/LOG_THERMAL_VERIFICATION.md) now corrects that
local axis defect in a separate component. Its checkpoint passed **361 tests**, with 18 audited
rate cases and nine full-transit prescribed-current transport cases. The matched
near-axis entropy-rate error drops from about 39% in the old representation to
8.85e-9 on the finest grid. Physical entropy and pressure work share the actual
conservative density and physical integration metric. Uniform-temperature
false magnetic response and full-electric-field precision remain measured
limitations. Bulk motion and the full conservative energy coupling are next;
these components have not yet replaced the nonlinear moving solver.

The [bulk thermal study](docs/BULK_THERMAL_VERIFICATION.md) now advances actual
changing density and both species entropies through one conservative particle
flux. That checkpoint passed **370 tests**, with 14 audited rates and ten audited compression and
full-transit trajectories. Finest translation density/temperature errors fall
below 1e-8; compression errors are around 1e-6. An initial face instability is
preserved and corrected. Conservative mechanical pressure work, full moving
magnetic coupling and the original heated-edge/formation requirements remain
unfinished; the prescribed-flow tests do not establish a complete plasma model.

The new [self-consistent fluid integration](docs/LOG_FLUID_NUMERICS.md) advances
momentum and both physical temperatures, with **376 passing tests**. Longer
tests caught a radial acoustic instability; the preserved diagnostic figure
shows the failure and a revised work pairing that completes a full contact
transit with a smaller step. The [radial-work analysis](docs/RADIAL_WORK_NUMERICS.md)
also exposes a remaining near-axis radial-force error, so this is still an
experimental component. The interrupted old remote comparison retains twelve
complete cases, one numerical failure and one missing fine contact.
[View the scientific diagnostics](outputs/log_fluid_diagnostics/report.png).

The new [common radial momentum study](docs/FACE_RADIAL_VERIFICATION.md) fixes
that near-axis force mismatch in a separate prototype. **382 tests pass**, nine
two-dimensional rate cases and five completed trajectories are independently
audited, and radial sound frequencies agree with the analytic cylindrical
reference. Remaining physical energy-transport error is measured explicitly;
finer checkpointed controls and full moving magnetic coupling remain ongoing.
[View the new convergence figure](outputs/face_radial_diagnostics/report.png).

The [joint entropy/momentum transport study](docs/TWO_ENTROPY_VERIFICATION.md)
now tests physical energy conservation through one shared flux. **385 tests
pass**, with 12 audited original-equation rates and 13 audited axial trajectories.
Both species remain accurate under refinement; remaining energy drift decreases
with the time step. The physical energy theorem applies to the declared nodal
integration, checked against denser continuous sampling. Cylindrical geometry
and compatible magnetic/current coupling still require integration.

[Explore the interactive 3D field showcase](outputs/frc_showcase/index.html)
or view its [high-resolution poster](outputs/frc_showcase/poster.png). It shows
archived current and magnetic fields from an initialized FRC verification case;
[the visualization notes](docs/SHOWCASE_3D.md) explain its data and limits.

![Waveform comparison](outputs/comparison/comparison.png)

Start with [the measured findings](docs/RESULTS.md),
[verification and sensitivity](docs/VERIFICATION.md),
[governing equations and assumptions](docs/PHYSICS.md), and
[architecture / stack decision](docs/ARCHITECTURE.md).

The continued development is tracked in [the active goal](docs/GOAL_STATUS.md).
The [external-engine comparison](docs/ENGINE_COMPARISON.md) records 17 audited
WarpX GPU fixtures and a Gkeyll source review. The existing research path remains
the core; the isolated WarpX installation is available for independent and
future kinetic comparisons, with its measured limits retained.
See [dynamic MHD verification](docs/MHD_VERIFICATION.md),
[curved closed-flux equilibrium checks](docs/EQUILIBRIUM_VERIFICATION.md), and the
[equations and boundary contract](docs/MHD_NUMERICS.md).
The [low-pressure verification](docs/POSITIVITY_VERIFICATION.md) distinguishes
positive reconstruction and conserved budgets from adequate flow resolution.
The [dynamic coupling study](docs/COUPLED_VERIFICATION.md) measures independent
force/advection/diffusion errors, external-source work and finite-domain error.
Its [energy and boundary contract](docs/COUPLED_NUMERICS.md) defines the scope.
The [open-boundary study](docs/OPEN_BOUNDARY_VERIFICATION.md) now verifies a
decaying magnetic exterior to infinity, including its energy beyond the grid.
The [strong magnetic-shock study](docs/SHOCK_VERIFICATION.md) exposes false field
reversal in the centered split method and measures an explicit numerical
diffusion option. Its heat is reported separately from physical Ohmic heating.
The [conservative coupling study](docs/CONSERVATIVE_VERIFICATION.md) advances
fluid and magnetic shock fluxes consistently through the same verified vacuum.
It measures a complete-update flux limiter, numerical interface heating and
whole-space energy error separately from exact fluid-cell conservation.
The [local-limiter study](docs/LOCAL_LIMITER_VERIFICATION.md) confines fallback
to troubled boundary stencils and restores interior shock accuracy. It retains
the newly exposed whole-space energy-error limitation explicitly.
The [compatible-energy study](docs/COMPATIBLE_ENERGY_VERIFICATION.md) tests
staggered magnetic energy, natural vacuum boundary work and optional RK3.
Its independent step audit separates spatial work consistency from remaining
time-integration error; these controls are explicit during verification.
The [boundary accuracy study](docs/BOUNDARY_VERIFICATION.md) adds an independent
cylindrical skin solution, full-field coil-pulse refinement and separate gas-end
and magnetic-matching controls. It measures an explicit first-order boundary
strip's accuracy and dissipation alongside the second-order comparison.
An [analytic spherical FRC reference](docs/FREE_BOUNDARY_REFERENCE.md) now tests
force balance, magnetic-vacuum matching and topology diagnostics. It is prescribed
equilibrium for free-boundary verification, not a formation result.
The [moving-material study](docs/MOVING_PLASMA_VERIFICATION.md) measures actual
free expansion, equilibrium drift, magnetic virtual force, source work and
interior flux interpolation error. Its [method](docs/MOVING_PLASMA.md) makes
the ideal component and mesh/axis approximations explicit.
The [resistive moving study](docs/RESISTIVE_MOVING_VERIFICATION.md) adds independent
spherical-decay, current, local-heat and stiff-time checks. Its
[method](docs/RESISTIVE_MOVING_NUMERICS.md) defines adaptive diffusion control and
atomic coupling to motion. The finest driven-flow comparison still differs by
63.8% between meshes despite much smaller energy and time-integration errors.
The [surface-resolution study](docs/SURFACE_VERIFICATION.md) adds independent
graded-mesh references, local force/thermal profiles and common physical-volume
comparisons. Its integrated velocity differences remain large; additional fine
time controls change velocity by 5.40%, then 2.76%. These use a different norm
from the historical 63.8% result. The [surface compatibility audit](docs/SURFACE_COMPATIBILITY.md)
derives a conflict between initial edge heating and a smooth zero-pressure
surface under the stated continuous-field assumptions. This restricts the
interpretation of these stress fixtures and does not prove that every weak
solution is impossible.
The [current-inertia study](docs/CURRENT_VERIFICATION.md) adds 35 complete
component cases and [physical regime checks](docs/TRANSPORT_REGIME.md). The
original argon reference has an ion inertial length about four times its plasma
radius, so its parameters do not justify ignoring Hall physics. A coarse
collisionless screening mesh also produces false reversal despite a small
energy residual; this artifact and the remaining fine-grid error are explicit.

The [linear Hall study](docs/LINEAR_HALL_VERIFICATION.md) adds 49 cylindrical
runs and 16 parallel-wave references, including both polarizations, finite
electron mass, damping and mesh/time refinement. This component verifies the
linear work pairing; its finite-mass extension to nonlinear moving flow remains
incomplete. The full suite now passes 293 tests.
The [nonlinear toroidal study](docs/TOROIDAL_MHD_VERIFICATION.md) adds
compressible rotation, toroidal induction and conservative angular momentum,
with 33 complete equilibrium/transport/wave/diffusion cases. Independent physical
field errors and heat convergence are measured. Physical free-edge coupling
remains unfinished.

The [Hall spatial study](docs/HALL_WORK_VERIFICATION.md) adds a verified
finite-pressure wave reference and 42 nonlinear/current/work audits. The
experimental operator conserves global Hall work but retains a 21.46% error
in the first annulus's induction rate. That earlier operator and its failed
data remain available as a historical control.

The [compatible Hall study](docs/COMPATIBLE_HALL_VERIFICATION.md) corrects that
axis error and adds actual nonlinear compressible evolution. All 63 spatial,
wave, nonlinear, time and resistive controls complete, with independent archive
energy replay. The first-annulus error reaches 0.00830% at nr128. Physical
magnetic reconstruction remains first order; fine wave errors are about 4%,
despite much smaller toroidal-coefficient and phase errors. The
[method](docs/COMPATIBLE_HALL_NUMERICS.md) states boundaries and reductions;
[separate species thermal transport](docs/TWO_TEMPERATURE_PLAN.md) is next.

The [electron-pressure study](docs/ELECTRON_PRESSURE_VERIFICATION.md) adds a
conservative current/heat-transport component with paired pressure induction,
electron energy and discrete physical entropy. Its 46 declared checks include
actual heat translation and nonlinear pressure-driven magnetic-field creation;
one deliberately strict solver failure is retained. The finest continuous
heat-transport error is 0.2672% after a full transit. This fixed-density component
still needs coupling to moving ion mass/energy and a unified boundary model.

The [moving-species study](docs/WARM_ADIABATIC_VERIFICATION.md) now verifies the
ideal common-velocity part of that coupling: shared conservative mass flow,
compression and separate temperatures. All 36 cases complete; finest continuous
electron-pressure errors are 0.561% in translation and 0.0593% in compression.
A stronger-case physical entropy decrease shrinks with refinement but remains
explicit. This ideal reference remains a separate limit of the new coupling.

The [warm Hall integration study](docs/WARM_HALL_VERIFICATION.md) adds 58 main
measurements and five successful supplemental controls, with conservative
electron-current enthalpy transport and a consistent zero-Btheta/free-Ez wall.
Finest fast-wave electron-pressure error is 0.347%, but the slow branch still
has 16.4% electron-pressure and 83.0% ion-pressure error. Three fine-grid
[pressure-split rejections](docs/WARM_PRESSURE_SPLIT_LIMITATION.md) preserve
negative prospective ion energy despite conserved global energy. An independent
continuum symbol also restricts that isolated split. Coupled nonideal integration,
axis accuracy and the full physical formation requirements remain active.

## Use the existing environment

```sh
.venv/bin/plasma-bench run --config configs/two_coils.json --stage 0 --out outputs/my_vacuum --viewer
.venv/bin/plasma-bench run --config configs/two_coils.json --stage 1 --out outputs/my_diffusion --viewer
.venv/bin/plasma-bench compare --peak-field 0.015 --out outputs/my_comparison
.venv/bin/plasma-bench sweep --out outputs/my_sweep
.venv/bin/plasma-bench view outputs/comparison/burst --gif
.venv/bin/python -m pytest -q
.venv/bin/python -m examples.hall_wave_measurements
.venv/bin/python -m examples.hall_parallel_measurements
.venv/bin/python -m examples.verify_hall
.venv/bin/python -m examples.toroidal_measurements
.venv/bin/python -m examples.verify_toroidal
.venv/bin/python -m examples.hall_work_measurements
.venv/bin/python -m examples.verify_hall_work
.venv/bin/python -m examples.compatible_hall_measurements
.venv/bin/python -m examples.compatible_hall_controls
.venv/bin/python -m examples.verify_compatible_hall
.venv/bin/python -m examples.electron_measurements
.venv/bin/python -m examples.verify_electron_pressure
.venv/bin/python -m examples.warm_measurements
.venv/bin/python -m examples.warm_controls
.venv/bin/python -m examples.verify_warm_adiabatic
.venv/bin/python -m examples.warm_hall_measurements
.venv/bin/python -m examples.verify_warm_hall
.venv/bin/python -m examples.verify_mhd
.venv/bin/python -m examples.verify_equilibrium
.venv/bin/python -m examples.verify_positivity
.venv/bin/python -m examples.verify_coupled
.venv/bin/python -m examples.verify_open_boundary
.venv/bin/python -m examples.verify_shocks
.venv/bin/python -m examples.verify_conservative
.venv/bin/python -m examples.verify_local_limiter
.venv/bin/python -m examples.verify_compatible_energy --reuse
.venv/bin/python -m examples.verify_spherical_frc
.venv/bin/python -m examples.verify_boundary
.venv/bin/python -m examples.moving_static_measurements
.venv/bin/python -m examples.moving_measurements --reuse
.venv/bin/python -m examples.verify_moving
.venv/bin/python -m examples.resistive_decay_measurements
.venv/bin/python -m examples.resistive_dynamic_measurements --reuse
.venv/bin/python -m examples.verify_resistive_moving
.venv/bin/python -m examples.surface_measurements --references
.venv/bin/python -m examples.surface_measurements --reuse
.venv/bin/python -m examples.surface_time_control
.venv/bin/python -m examples.verify_surface
.venv/bin/python -m examples.current_measurements --reuse
.venv/bin/python -m examples.transport_study
.venv/bin/python -m examples.verify_current
```

The boundary verifier assembles completed boundary-study records. To regenerate its
pulse, skin and end measurements, run the `examples.boundary_pulse_measurements`,
`examples.cylindrical_skin_measurements` and `examples.boundary_end_controls`
modules. The separately recorded fine shock uses
`examples.boundary_order_probe.run(384, .2, "outputs/boundary_verification/fine_shock")`.
These fine-grid studies can take many minutes; partial records are not results.
The moving measurement/assembly commands validate exact source/configuration and
lossless field hashes before reusing completed cases. The resistive commands
reproduce the uniform extension; the surface commands reproduce the graded
study; the final three reproduce the current-inertia and transport assessment.
Historical studies retain their recorded
source fingerprints; current production changes do not relabel old outputs.
The finest material meshes also take several minutes.

Open [vacuum playback](outputs/vacuum/viewer.html) or
[conducting-plasma playback](outputs/comparison/burst/viewer.html) in a browser.
The standalone viewers include saved frames, play/pause, frame stepping and a
slider; no server or internet connection is required. This is playback of real
computed fields. Edit configuration and rerun to change the experiment.

## Install on another machine

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt -e .
.venv/bin/python -m pytest -q
.venv/bin/python -m examples.hall_wave_measurements
.venv/bin/python -m examples.hall_parallel_measurements
.venv/bin/python -m examples.verify_hall
.venv/bin/python -m examples.toroidal_measurements
.venv/bin/python -m examples.verify_toroidal
```

Python 3.11+ is supported by the package; the exact delivered environment is
Python 3.13 on macOS ARM64. The lock lists the versions actually tested here;
installing `-e '.[dev]'` instead uses the declared compatible version ranges.
CUDA and a GPU are not required. Source fields and sparse factorizations are
reused. An RTX backend is a future profiling-driven option.

## Script an experiment

```python
from dataclasses import replace
from plasma_bench.config import Experiment
from plasma_bench.coils import Coil
from plasma_bench.waveforms import DampedSine
from plasma_bench.simulation import run

base = Experiment()
wave = DampedSine(amplitude=1000, frequency=50_000, decay=20e-6)
config = replace(base, coils=(
    Coil(radius=0.061, z=-0.055, waveform=wave),
    Coil(radius=0.061, z=0.055, waveform=wave),
))
result = run(config)
result.save("outputs/custom")
print(result.metadata["summary"])
```

All inputs are SI; temperatures are explicitly named in eV. Coils support
radius, axial position, effective turns and polarity. Waveforms include DC,
Gaussian, cosine-ramped unipolar, damped sinusoid, smooth finite sinusoidal burst
and sampled data. Analytic field evaluation also accepts script callbacks;
sample them before persistent simulations so output remains reproducible.

`plasma-bench init configs/new.json` writes the full configuration. Material
species mass, charge, density, temperatures, ionization fraction and prescribed
conductivity are explicit. Gas pressure is **not** converted into density or
conductivity without an ionization/transport closure. See the research queue.

## Repository map

- `geometry.py`, `config.py`: typed chamber, outer domain, plasma and experiment.
- `coils.py`, `waveforms.py`, `electromagnetics.py`: solver-independent sources.
- `numerics.py`, `plasma.py`: compatible flux operators, mixed diffusion/vacuum solver.
- `mhd*.py`: developmental conservative dynamics and independent analytic benchmarks.
- `coupled*.py`, `coupling_mesh.py`: gas/magnetic splitting, prescribed coils and genuine exterior vacuum.
- `open_boundary.py`: decaying magnetic multipoles and full-space magnetic energy.
- `shock_reference.py`, `shock_benchmarks.py`: independent cylindrical magnetic shocks and actual-kernel comparisons.
- `magnetic_dissipation.py`: optional declared compressive numerical diffusion, with separately accounted heat.
- `conservative_coupled.py`, `vacuum_flux.py`, `vacuum_trace.py`: unsplit conservative fluid/CT evolution, interface traction and magnetic vacuum.
- `flux_update.py`, `flux_limiting.py`: shared flux assembly, local complete-update fallback and explicit global comparison.
- `magnetic_quadrature.py`, `vacuum_energy.py`: positive staggered energy and natural magnetic-vacuum work.
- `moving_geometry.py`, `moving_field.py`: evolving material geometry, magnetic virtual work and exact infinite exterior.
- `moving_advection.py`, `moving_plasma.py`: ideal material flux labels, adiabatic gas and actual free-surface motion.
- `material_resistivity.py`, `resistive_substeps.py`, `resistive_moving.py`: consistent current, physical diffusion/heat, adaptive accuracy control and atomic resistive/material composition.
- `current_inertia.py`: separate fixed-material linear current response, physical electron inertia, drift kinetic energy and true exterior vacuum.
- `coulomb_transport.py`: conditional classical collision coefficients and explicit local carrier/Hall/inertia/thermal scale checks.
- `diagnostics.py`, `topology.py`: validity scales, energy, O/X points, closed flux.
- `simulation.py`, `results.py`: orchestration and reproducible numeric archives.
- `experiments.py`, `visualization.py`, `cli.py`: sweeps, plots and playback.
- `tests/`, `examples/`, `docs/`: verification, runnable experiments and research notes.

Every Python module has an adjacent Markdown contract. Each saved run contains
`fields.npz` (numeric arrays, no pickle) and `metadata.json` (full configuration,
units, source fingerprint, versions, warnings, topology and summary).
Scalar histories resolve every integration step; field/topology snapshots use
`output_stride`, and playback may subsample these again for file size.

Regenerate the deliverables with `python -m examples.build_demo`; generate the
convergence report with `python -m examples.verify`. These run after the test
suite. The 27-case phase map is exploratory, not a validated formation boundary.

## Next physics milestone

Second-order reconstruction reduces measured numerical diffusion and axis-force
error; local admissibility limiting resolves the identified negative face-pressure
failures without cell floors. A fixed gas-filled chamber now couples moving fluid
to prescribed coils and genuine surrounding magnetic vacuum with a measured
source-work budget. A verified multipole boundary extends the induced magnetic
field to infinity and counts its exterior energy, without adding fluid there.
An unsplit conservative vacuum-coupled reference now preserves the tested
perpendicular shock's frozen-flux relation while measuring interface and global
energy errors. Local complete-update fallback preserves distant interior
accuracy, and compatible spatial work plus a quadratic RK audit separate the
remaining energy mechanisms. A fixed boundary strip improves shock energy but
has a measured driven-flow accuracy cost. The new ideal material reference
provides actual free-surface evolution and independent expansion/force checks.
Its resistive extension verifies scalar diffusion and heating, but resolved
driven flow, a defensible heated-edge model, spatial coil coupling, realistic
ends and transport/Hall validity remain required. Surface grading preserves
the uniform default and makes refinement explicit; it does not by itself
resolve the flow or repair the conditional thermal incompatibility.
The finite-inertia current component is verified separately and is not yet
coupled to moving gas. The original argon scales motivate Hall/toroidal and
electron/ion thermal extensions; [the next integration contract](docs/EXTENDED_PLASMA_PLAN.md)
retains their energy, vacuum and model-validity requirements.
The coil-driven Stage 2 entry point still raises an explicit `NotImplementedError`
until the physical formation requirements are verified. The dynamic references
are separately scriptable through `plasma_bench.mhd.AxisymmetricMHD` and
`plasma_bench.coupled.VacuumCoupledMHD`; the latter's gas energy excludes magnetic
energy and must not be interchanged with the former's total-MHD state.
Cold-rarefaction pressure accuracy and curved-equilibrium drift remain measured
numerical limitations to address alongside that integration.

To run the conservative numerical coil fixture with the infinite magnetic exterior:

```python
from plasma_bench.coupled_benchmarks import coil_pulse
from plasma_bench.conservative_benchmarks import from_split_initial
model = from_split_initial(coil_pulse(16, vacuum_boundary="multipole", multipole_modes=12))
model.advance(0.1)
print(model.diagnostics())
```

These unit-scale SI verification parameters are not a physical argon experiment.
The matching surface must enclose the entire fluid source sphere; mode and mesh
refinement remain required. `vacuum_boundary="fixed"` retains the historical
finite-wall comparison. Neither option yet provides an evolving free gas surface.
The original split reference is available directly as `coil_pulse(...)`.
`ConservativeState` includes fluid-cell magnetic energy, whereas `CoupledState`
contains gas energy only; `from_split_initial` makes that conversion explicitly
for unadvanced initial fixtures. The conservative solver's whole-space energy
error is diagnosed, not forced to zero. See its [numerical contract](docs/CONSERVATIVE_COUPLING.md).
The default `flux_limiter="local"` replaces shared fluxes around inadmissible
cells and rechecks neighbors. `flux_limiter="global"` retains the earlier
comparison, and `positivity_flux_limiting=False` retains unprotected failures.
Neither mode repairs cell energy or guarantees every low-order proposal is valid.

The compatible-energy comparison is selected when constructing the initial
fixture, before any state is advanced:

```python
model = from_split_initial(
    coil_pulse(16, vacuum_boundary="multipole"),
    energy_quadrature="staggered", boundary_poynting="variational", time_order=3,
)
model.advance(0.1)
```

Do not change an evolved state's magnetic-energy convention. All three selectors
retain historical choices for controlled comparison. Neither a smaller energy
residual nor a higher time order establishes a converged free-plasma experiment.

To reproduce the initialized strong-shock verification fixture:

```python
from plasma_bench.shock_benchmarks import standing_shock
model, reference = standing_shock(48, shock_diffusion=1.0)
model.advance(0.02)
print(model.diagnostics())
```

`shock_diffusion` is numerical stabilization, off by default, not material
resistivity. The centered `0.0` comparison can generate false magnetic reversal;
neither option should be used to assert FRC formation. The independent exact
shock is valid only in the stated interior annulus and comparison time.
`midpoint_method="picard"` caches the vacuum solve and is now the faster default;
`"direct"` retains the original solve for comparison of the same equations.
Both require converged nonlinear residuals and leave accepted state intact on
failure. Read the [shock equations and limitations](docs/SHOCK_NUMERICS.md).

To evolve the independently known freely expanding gas sphere:

```python
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.moving_plasma import MovingPlasma

domain = MovingDomain(8, 16, radius=1.0, matching_radius=2.0)
model = MovingPlasma(domain, initial_a=0.0,
                     pressure=lambda r, z: 1.0-r*r-z*z,
                     max_timestep=0.001)
model.advance(0.1)
print(model.diagnostics())
```

Here `initial_a=psi/r²` includes the regular axis limit and is measured in tesla.
The model has no exterior gas, resistivity, shock capture or remapping. A moving
material mesh can distort; preserved vertex flux does not prove exact frozen
flux everywhere between vertices. `examples.moving_cases.MovingCase` supplies
initialized FRC and uniform-bias fixtures with fully recorded mathematical inputs.
They do not establish formation, species validity or confinement.

The moving resistive fixture is separately scriptable:

```python
from examples.resistive_records import ResistiveCase

model = ResistiveCase("example", "source", radial_cells=8).build()
model.advance(0.04)
print(model.diagnostics())
```

This initializes a known reversed-field equilibrium and applies a uniform bias
ramp. It tests resistive relaxation and source coupling, not formation from an
unreversed plasma. `ResistiveMovingPlasma` accepts explicit scalar resistivity in
ohm m; its diffusion tolerance controls local field/current/heat estimates, not
global physical accuracy. The [complete study](docs/RESISTIVE_MOVING_VERIFICATION.md)
retains the large moving-flow mesh dependence and source/time/boundary controls.
`ResistiveCase(..., radial_grading=2.0)` clusters initial material and magnetic
vacuum shells at the shared surface. The default is `1.0`; the new uniform
runs reproduce the historical saved state arrays exactly. See the
[surface method](docs/SURFACE_NUMERICS.md) for its measured accuracy and limits.

Stage 0/1 run metadata now includes conditional Coulomb estimates from the
supplied species and temperatures. Their reference current is initially zero;
the local `plasma_bench.coulomb_transport.regime` function accepts actual current,
physical structure length and variation time for a separate assessment.
These estimates never overwrite prescribed conductivity or establish validity
merely because a warning threshold is not crossed.

[Open physics and numerical questions](docs/RESEARCH_QUESTIONS.md) include the
outer vacuum boundary, transport, pre-ionization, neutral effects, objective FRC
criteria and synthetic probe responses. All hardware remains prescribed coil
current; the project contains no pulsed-power construction design.
