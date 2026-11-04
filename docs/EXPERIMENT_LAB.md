# Local experiment lab

Fieldwork is the first unified local app for this research workspace. A shared
experiment registry, run library, launch/status controls and replay shell can
support future experiments while keeping each solver's equations and validity
separate. It does not require a universal plasma solver or a cloud service.

## Use

```sh
.venv/bin/python -m pip install -e '.[electrostatic]'
.venv/bin/python -m examples.lab_app
```

Open `http://localhost:8792`. The server binds loopback only. The optional
electrostatic dependency is Gmsh 4.15.2; the Poisson/particle path uses CPU64.
Existing FRC server/jobs are independent and do not need to be stopped.

- **Positive cage:** configure physical initial plasma, electron source and
  voltage drive; launch a full 3D electron/ion Vlasov–Poisson experiment. Saved
  frames show actual particles, potential cut, energy ledger and input values.
  Orbit/zoom, timeline and display-count slider affect visualization only.
- **FRC formation:** existing 2D and 3D saved runs are accessible through their
  original scientific viewer, including current volumes, field lines and ions.
  This first adapter does not launch/adopt FRC campaigns or refresh their exports.
- The app owns at most one IEC child. Stop writes a request to that child only;
  partial scientific data stay. Stopping the web server does not stop simulations.

The initial launcher supports 1–2,000 ns and at most 200,000 timesteps, with
explicit particle/mesh bounds. These are resource limits, not a physical validity
window. A 50 ns run follows electron startup, not a full POPS cycle. Initial
plasma density is specified directly; no pressure-to-plasma conversion is implied.

## Structure

`plasma_bench/lab_catalog.py` contains experiment definitions, physical input
schema/validation and the FRC read adapter. `examples/lab_app.py` serves the app
and maps validated IEC inputs to a fixed solver command. New model adapters can
reuse navigation, run manifests and viewing without inheriting FRC quasi-neutral
electrons or electrostatic assumptions. There is no arbitrary shell/path API.

`visualizations/lab/` contains buildless HTML/CSS/JS, with no external scripts.
Run URLs retain `?run=ID`. Data remain in `outputs/`; no outputs are committed.
Each launched IEC run retains source copies/hashes, mesh, input configuration,
atomic progress/history, actual-time frame JSON, full 3D state NPZ and final state.
No displayed particle sample substitutes for the full saved physical population.

## End-to-end checks (2026-10-21)

An actual 20 ns static-bias simulation was launched and completed from the form.
Browser replay, time scrubbing and zero/reduced dot counts were checked: physical
counts remained unchanged. A 200 ns owned run was stopped during preparation,
retaining its initial state and cancelled status. Both existing 2D formation and
the complete 0–40 µs 3D hybrid viewer were loaded through the adapter, including
native half-float current volumes. Preparing a new run clears old-run graphics.
The final source-driven replay is recorded in [IEC_3D.md](IEC_3D.md).

This is a first usable research app, not completed common run control for every
solver. FRC campaigns keep their existing source-frozen tools. Scientific limits
remain in each model's scope; validation details for the new solver are in
[IEC_3D.md](IEC_3D.md).
