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
