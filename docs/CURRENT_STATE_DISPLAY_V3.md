# Updated magnetic and signed-current display

Open [the interactive viewer](../outputs/current_state_display_v3/index.html) or
[the inspected scientific PNG](../outputs/current_state_display_v3/current_state_3d.png).
The optional localhost viewer is at http://127.0.0.1:8788. It was started in
session 60135; it is a presentation server, not a scientific worker.

The latest evolved state remains **step 17 at 28.3333 µs**. The subsequent
entropy, spatial-refinement and coupled-rate diagnostics have not advanced the
trajectory. The filled periodic verification cylinder is 70 mm in diameter and
240 mm long with a 5 mT axial bias. These are verification dimensions, not a
recommended laboratory design. Step 18 rejected; 100 µs was not reached.

The viewer adds signed radial, azimuthal and axial current selection alongside
total magnitude, with symmetric signed legends and shared or per-frame scales.
Four actual saved times, linked rotation/zoom, camera presets, cutaway, opacity
and PNG export remain. Magnetic lines use the complete field including bias;
their deflection and the physical aspect ratio are not exaggerated.

The scientific image adds meridional axial-current and `Bz − 5 mT` plots. Their
colors are signed physical values, unlike the illustrative optical blending in
the interactive volume. Full B still determines every magnetic line.

At the final saved time, display-sampled peak |J| is **835.178 kA/m²**, signed
Jz spans **−835.178 to +500.698 kA/m²**, and minimum Bz is **+4.92882 mT**.
The axial field remains positive; no field reversal or FRC formation is shown.
These extrema are display samples, not certified extrema of the physical field.

## End-to-end verification

Both complete display builds finished with collected exit code 0. HTML, PNG,
packing metrics and manifests match byte for byte in fresh processes. All 560
recorded inputs were rechecked, including the preceding 200-source preservation
record and completed newer scientific replay evidence. Embedded full magnetic
paths equal the float64 archive exactly. Every signed texture value equals its
archived component cast to float32, and the original magnitude texture remains
unchanged. Maximum conversion error divided by component peak is 5.333e-8.

The 2880 × 2070 PNG was visually inspected. Assembled JavaScript passes Node
syntax checking. Browser interaction remains **unverified** because the Mac is
locked; the app queued the browser and image tabs. No unit tests or solver runs
were performed for this update. The user's earlier poster is unchanged.

See [the complete replay record](../outputs/current_state_display_v3_verification/summary.json)
and the [build manifest](../outputs/current_state_display_v3/manifest.json).
The HTML, PNG and manifest served over localhost also match disk exactly;
see [the delivery check](../outputs/current_state_display_v3_verification/served.json).
The inherited reconstruction, interpolation and tracing limits remain those in
[TRAJECTORY_DISPLAY.md](TRAJECTORY_DISPLAY.md). In particular, off-grid display
interpolation errors and the axisymmetric scope are not improved by this
presentation. The new slice plots use nearest-pixel cell-center samples and
half-pixel edge extension; no physics smoothing or numerical diffusion is added.

## Rebuild and continuation

```sh
PYTHONDONTWRITEBYTECODE=1 VECLIB_MAXIMUM_THREADS=1 \
MPLCONFIGDIR=/private/tmp/plasma-bench-mpl MPLBACKEND=Agg \
.venv/bin/python -m examples.current_state_display_v3 --output outputs/my_new_display
```

Completed destinations and used sources are protected from replacement. Add
`--serve` only if port 8788 is free, or choose another port. The HTML also works
offline. Read [the scientific continuation boundary](COUPLED_FLUID_REFINEMENT_PROGRESS.md)
before resuming simulation work. Full goal and acceptance gates remain open.
