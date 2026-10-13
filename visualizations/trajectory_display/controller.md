# Saved-state viewer controls and export

Synchronizes camera rotation, zoom, cutaway and saved-frame changes between two
WebGL views. Saved times are discrete and never interpolated. Auto orbit is off
initially; it moves only the camera. Current scaling defaults to the maximum
over all four snapshots and has an explicit per-frame inspection option.

The opacity slider affects optical presentation only. Every frame change
updates the sampled current, total B and minimum axial B values with SI units.
`poster` creates a labeled 2800×1900 PNG with the selected time, scales, physical
dimensions, axisymmetric scope and trajectory rejection. Local server export
posts only that PNG to `/save-poster`; standalone HTML downloads the same image.

UI verification uses actual browser interaction and exported image inspection,
without unit tests. Scientific data reconstruction/replay is separate.
