# Virtual plasma bench

An axisymmetric plasma simulator (Python, NumPy/SciPy) used to design a real
benchtop field-reversed configuration (FRC) experiment: a field-reversed theta
pinch in a ~15 cm × 1 m Pyrex tube with H2/D2 fill.

- **Current status and plan:** [docs/GOAL_STATUS.md](docs/GOAL_STATUS.md)
- **First-cut fields, timing and bank numbers:** [docs/FRC_DESIGN_WINDOW.md](docs/FRC_DESIGN_WINDOW.md)
- **What exists and its limits:** [docs/PROJECT_OVERVIEW.md](docs/PROJECT_OVERVIEW.md)
- **Working rules:** [AGENTS.md](AGENTS.md); issues are tracked with `bd` (epic `frc-ibl`)

This is research software. No result here has yet been benchmarked against a
formation experiment, and nothing yet demonstrates FRC formation in this geometry.

## Install and run

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt -e .
.venv/bin/python -m examples.frc_design_window          # design window, seconds
.venv/bin/python -m examples.material_experiment configs/material_tube_bias_pulse_20mt.json outputs/my_run
```

Python 3.11+; the delivered environment is Python 3.13 on macOS ARM64. No GPU
is required. `plasma-bench init/run/compare/sweep/view` exposes the original
stage 0–1 vacuum and resistive-diffusion workflows.

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
- `frc_scaling.py`: first-cut field-reversed theta-pinch design relations.
- `tests/`, `examples/`, `docs/`: historical tests, runnable experiments and research notes.

Every Python module has an adjacent Markdown contract. Each saved run contains
`fields.npz` (numeric arrays, no pickle) and `metadata.json` (full configuration,
units, source fingerprint, versions, warnings, topology and summary).
Scalar histories resolve every integration step; field/topology snapshots use
`output_stride`, and playback may subsample these again for file size.

Simulation outputs are not tracked in git; each example regenerates its own
directory under `outputs/`.
