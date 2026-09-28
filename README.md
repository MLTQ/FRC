# Virtual plasma bench

A scriptable, verified axisymmetric research sandbox. It currently implements
analytic vacuum electromagnetics and stationary resistive plasma induction with
a coupled nonconducting exterior. A developmental conservative MHD kernel now
simulates compression, axial motion, pressure/energy evolution and resistive
heating in explicit verification domains. The coil-driven formation workflow
now has a separate fixed-chamber coupling reference with moving gas, prescribed
coils and an elliptic magnetic vacuum. Free plasma interfaces, appropriate end
conditions and physical transport remain unfinished. Stability and particle
confinement have not been established.

![Waveform comparison](outputs/comparison/comparison.png)

Start with [the measured findings](docs/RESULTS.md),
[verification and sensitivity](docs/VERIFICATION.md),
[governing equations and assumptions](docs/PHYSICS.md), and
[architecture / stack decision](docs/ARCHITECTURE.md).

The continued development is tracked in [the active goal](docs/GOAL_STATUS.md).
See [dynamic MHD verification](docs/MHD_VERIFICATION.md),
[curved closed-flux equilibrium checks](docs/EQUILIBRIUM_VERIFICATION.md), and the
[equations and boundary contract](docs/MHD_NUMERICS.md).
The [low-pressure verification](docs/POSITIVITY_VERIFICATION.md) distinguishes
positive reconstruction and conserved budgets from adequate flow resolution.
The [dynamic coupling study](docs/COUPLED_VERIFICATION.md) measures independent
force/advection/diffusion errors, external-source work and finite-domain error.
Its [energy and boundary contract](docs/COUPLED_NUMERICS.md) defines the scope.

## Use the existing environment

```sh
.venv/bin/plasma-bench run --config configs/two_coils.json --stage 0 --out outputs/my_vacuum --viewer
.venv/bin/plasma-bench run --config configs/two_coils.json --stage 1 --out outputs/my_diffusion --viewer
.venv/bin/plasma-bench compare --peak-field 0.015 --out outputs/my_comparison
.venv/bin/plasma-bench sweep --out outputs/my_sweep
.venv/bin/plasma-bench view outputs/comparison/burst --gif
.venv/bin/python -m pytest -q
.venv/bin/python -m examples.verify_mhd
.venv/bin/python -m examples.verify_equilibrium
.venv/bin/python -m examples.verify_positivity
.venv/bin/python -m examples.verify_coupled
```

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
source-work budget. The exterior contains no fluid atmosphere. Next, address
free interfaces, realistic end conditions, source-driven skin-layer resolution
and magnetic-shock fidelity, alongside conductivity/Hall applicability.
The coil-driven Stage 2 entry point still raises an explicit `NotImplementedError`
until the physical formation requirements are verified. The dynamic references
are separately scriptable through `plasma_bench.mhd.AxisymmetricMHD` and
`plasma_bench.coupled.VacuumCoupledMHD`; the latter's gas energy excludes magnetic
energy and must not be interchanged with the former's total-MHD state.
Cold-rarefaction pressure accuracy and curved-equilibrium drift remain measured
numerical limitations to address alongside that integration.

[Open physics and numerical questions](docs/RESEARCH_QUESTIONS.md) include the
outer vacuum boundary, transport, pre-ionization, neutral effects, objective FRC
criteria and synthetic probe responses. All hardware remains prescribed coil
current; the project contains no pulsed-power construction design.
