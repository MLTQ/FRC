# formation_display.py

Interactive 3D pages for formation runs from `examples/nucte_formation.py`
(coupled drive), using the saved-field viewer in
`visualizations/material_fields`.

```bash
python -m examples.formation_display OUT RUN[::Label] ...
```

It writes four things:
- `OUT/viewer.html`: one 3D viewer for all runs.
  - A run selector switches between experiments, keeping the camera and the
    simulated time (nearest saved state).
  - `viewer.html?run=NAME` opens a given run.
- `OUT/<run>/manifest.json` and `OUT/<run>/frame_NNNN.json`: the run's viewer
  data, plus one file per saved state.
  - Frames load on demand (cached and prefetched) and play at up to 20
    frames/s.
  - Each is an actual saved state, never an interpolation. Use
    `snapshot_s=1e-7` in nucte_formation for smooth playback.
- `OUT/index.html`: trapped flux, separatrix length and midplane radius
  against time for up to six runs, plus a table linking into the viewer.
- `OUT/runs.json`.

Serve `OUT` over HTTP, for example with the `formation-display` entry in
`.claude/launch.json`, which serves `outputs/` on 127.0.0.1:8791. Then open
`/formation_display/index.html`.

## What is shown
Every frame is a saved 1 µs snapshot; nothing is evolved or interpolated in
time.

**Flux.** The run saves total flux on the gas nodes only. The builder rebuilds
the run's coupling mesh and section coil model from `setup.json`, then:
1. subtracts the coil flux of the saved section vacuum currents;
2. extends the induced flux into all vacuum (outside the wall, inside an inner
   tube, beyond the ends) with the run's own exterior solve,
   `CouplingMesh.harmonic_extension`;
3. adds the coil flux back.

The gas-node flux is reproduced to round-off (~1e-16), which is checked and
printed per run.

**Field lines** are flux contours in four meridional planes, rotated copies
for ×2 and ×4.
- Eleven open levels are equally spaced in flux, so line density follows
  field strength.
- Four closed levels lie between the O-point and the separatrix.
- Colour is |B| from the same flux. The colour scale is set by the field in the
  gas, which keeps line colours off the coil filaments.

**Current.** J_θ = −Δ*ψ/(μ0 r) from second differences, shown on gas nodes
only and rendered as a WebGL volume. J_r = J_z = 0, since the model has no
toroidal field.

**Separatrix and FRC numbers** come from `nucte_compare.summarize`
(`closed_region` or `annular_closed_region`). A run counts as formed while it
meets the design-scan criterion: closed, and shorter than 80% of the tube.

**Coil groups** are the main coil and each AuxCoil, each in its own colour.
Their per-frame labels give the actual coil current (vacuum + image) and its
vacuum Bz at the midplane tube wall.

## Overview page
The six series colours are the dataviz reference dark palette (blue, orange,
aqua, yellow, magenta, green). They are validated on the `#10161e` chart
surface: worst adjacent CVD ΔE is 8.4. Inner-tube runs are dashed as a
secondary encoding. There is one y-axis per chart, a hover crosshair and
tooltip, a legend that toggles series, and a table view.
