# Linked 3-D saved-field canvases

Projects the saved physical coordinates with one common rotation, perspective
and uniform length scale. Drag/zoom affect both canvases. Current colors support
magnitude (the default) and signed cylindrical components. Explicit per-state/shared scales
retain physical units and quantitative labels; transparency is display-only.

Playback and the slider select saved states without temporal interpolation.
The moving outline is supplied by the actual material mesh. Magnetic paths and
current samples come from the read-only builder; the renderer does not alter
scientific fields, evolve the simulation or infer topology/confinement.

No external packages or browser storage are used. Inline sparse displays remain
usable as one offline HTML file. Dense manifests fetch local JSON frame files.

Configured coil experiments may supply static chamber and coil paths, rendered
in grey and copper respectively, plus a per-state prescribed centre-field label
and an accurate status note. Missing optional fields retain the original compact
display behavior. These outlines are physical geometry, not adjusted aspect ratios.
Exactly zero field/current frames use the zero color without dividing by a zero
range; their displayed physical maxima remain exactly zero.
An optional initial_zoom fits the long tube more clearly; reset restores it.
Zoom always scales every coordinate equally.
An optional motion label shows the saved midplane and maximum boundary radius,
so a local pinch is distinguishable from expansion elsewhere. It reports actual
geometry supplied by the builder; the browser does not infer confinement.

Optional coil_groups replace the copper-only outlines with named colored groups.
Per-frame drive_groups show each group's vacuum-centre contribution alongside
the total. Text is inserted as textContent, with no HTML evaluation.

Optional solver_dt_s and solver_steps expose the actual integration step in ns,
separately from frame count and the gap to the preceding displayed state. Missing
metadata shows only the known frame count; no step size is inferred from animation.
Optional solver_completed_steps distinguishes a terminal partial run's evolved
step count from its requested count; complete runs retain the compact label.

Dense manifests may replace each frame's arrays with a relative `url`. The viewer
loads that file, caches at most five complete frames and prefetches one neighbor.
The visible and requested frames are protected from eviction during fast scrubs.
Scrubbing is briefly debounced; request identifiers keep stale loads from moving
the slider backward. A failed fetch is shown explicitly. No frame is synthesized.
All field paths/cloud points retain the builder's existing detail. Play waits for
each actual frame rather than skipping scientific states when drawing is slower
than the chosen1/5/10/20fps. Previous/next buttons move one displayed state.
`initial_frame` and `playback_fps` choose the starting state and playback rate;
`display_interval_ns` labels the sampling target separately from actual gaps.
The gap label refers to the preceding frame in time, including during reverse stepping.
Dense time labels use ns with two decimal places, exposing1.75ns-step rounding.
Field-line segments are batched into256 color levels per draw, avoiding tens of
thousands of canvas stroke calls. Every path vertex and actual stored field value
is retained; this quantizes only rendered color, not geometry or scientific data.

Detailed displays optionally attach float32 current textures to each frame.
These load with the original line/cloud data in the same bounded cache. Manifest
extrema include dense basis samples and take precedence over the old cloud-only
extrema. A WebGL volume layer shares an orthographic camera with the line canvas;
the original perspective is unchanged for older builders. All lengths retain
one common scale; a 50 mm screen-plane bar makes magnification explicit.
The detailed palette matches the volume shader. Signed volume colors are
blue/white/red. Geometry overlays are optional and fit/focus adjust only camera.
1/2/4 azimuthal copies give 96/192/384 paths for the current 96-path source;
this is exact rotational symmetry, not new field integration or new 3-D modes.
Opacity and threshold control only the illustrative volume transfer function.
The current cloud is the visible fallback if WebGL2 cannot initialize.

Detailed magnetic labels use mT, and large current legends use MA/m². Saved data
remain in tesla and A/m²; shared and per-state scales use the same unit conversion.

Detailed playback starts with a shared color scale, so the pulse brightens and
fades against fixed extrema; per-state scaling remains an explicit option.
The peak-current metric chooses A, kA or MA/m² according to its actual value.

## Inline volumes, separatrix and labels for formation runs
These optional fields are used by `examples/formation_display.py`. Missing
fields keep the earlier behaviour.

**`frame.volume.b64`**
- The current texture is carried inline as base64 little-endian float32, so a
  single HTML file works offline and over `file://`.
- With `channels` < 4, each texel holds that many values. `slots` places them in
  the RGBA texel (Jr, Jθ, Jz, support); the default `[1, 3]` is (Jθ, support).
  Other slots are zero.
- A size mismatch raises the same "Incomplete current volume" error as a fetched
  volume.

**`frame.separatrix`**
- Paths drawn in gold on both canvases, regardless of the geometry checkbox.
  An optional `#separatrix` checkbox, added by the formation builder and checked
  by default, hides or shows them.
- The colour can be set with `DATA.separatrix_color`.

**Text fields**
- `frame.summary`: a text line shown in `#frame-summary` when the page has
  that element.
- `DATA.drive_label` and `DATA.motion_label` replace the "Prescribed
  vacuum-centre field" and "Plasma radius — middle" prefixes.
- `DATA.temperature_digits` switches temperatures to that many significant
  digits. A null temperature shows "—".

## Several runs in one page
With `DATA.runs` = [{name, label, manifest}], `DATA.run` and a `#run` select,
the viewer switches runs in place:
1. It fetches the run's manifest.
2. It replaces DATA's fields and starts new frame caches. Loads still in flight
   finish into the old caches and cannot leak into the new run.
3. It moves to the saved state nearest the current simulated time and keeps
   the camera, colour-scale mode and current component.
4. Playback resumes if it was playing.

`#run-title`, `#run-description`, the document title and the scope note follow
`DATA.title`, `DATA.description` and `DATA.scope`. The URL becomes `?run=NAME`,
and loading with that parameter opens that run.

## Particle trails
With `DATA.particles` (a url) and a `#particles` checkbox, both added by
`examples/formation_display.py` when a run has `particles.json` from
`examples/particle_tracking.py`, the field canvas draws test-particle orbits.

**Data.**
- `{species, t0_us, sample_us, colors, sets: [{label, trails}]}`, where each
  trail is `{class, lost_us, surface, points}`.
- Points are viewer coordinates sampled every `sample_us` from `t0_us`. A lost
  trail ends with its exact loss point.

**Drawing.**
- For the current saved time t, each trail shows its last 2 µs (`TRAIL_US`),
  in its class colour, ending in a dot.
- A lost trail freezes at its loss point with a red ×, faded.
- Nothing is drawn before `t0_us`.
- Every vertex is an integrated orbit position. Nothing is interpolated or
  synthesized.

**Loading.**
- The file loads on first use and is replaced on a run switch.
- A run without orbits disables the checkbox.
- `#particle-set` selects the set (for example with and without the
  quadrupole) and is hidden when there is only one.
- `#particle-legend` names the classes, using textContent only.

The field lines remain those of the saved axisymmetric field. A coil field
superposed in the tracking (the quadrupole) is not drawn.

## Cartesian hybrid runs
`hybrid3d_display` supplies `no_copies`, fetched `volume3d` half-float textures,
and per-frame `ions` samples. The texture's `value_scale` restores A/m² during
rendering; legends always retain the physical extrema. Field-line copies are
hidden for these genuinely 3D fields. Ion dots are independent samples and
are never joined into trails. Model badge, temperature label and source notes
follow the selected manifest. Current components with nonzero run extrema are
enabled; unavailable components reset to magnitude on a run switch.
An optional ion-legend names the dot speed scale and saturation, separately
from the field-line |B| legend, and identifies each frame as a new sample.
