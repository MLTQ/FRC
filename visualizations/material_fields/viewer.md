# Linked 3-D saved-field canvases

Projects the saved physical coordinates with one common rotation, perspective
and uniform length scale. Drag/zoom affect both canvases. Current colors support
magnitude (the default) and signed cylindrical components. Explicit per-state/shared scales
retain physical units and quantitative labels; transparency is display-only.

Playback and the slider select saved states without temporal interpolation.
The moving outline is supplied by the actual material mesh. Magnetic paths and
current samples come from the read-only builder; the renderer does not alter
scientific fields, evolve the simulation or infer topology/confinement.

No external packages, network requests or browser storage are used. The single
offline HTML file remains usable independently of the optional local server.

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
