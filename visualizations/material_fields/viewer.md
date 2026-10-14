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
