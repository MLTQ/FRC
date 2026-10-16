# Axisymmetric current volume renderer

`createVolumeRenderer(canvas)` owns a WebGL2 texture/program and returns `draw`.
The shared viewer supplies one saved frame, current component, physical scale,
and the same orthographic camera used for the field lines. No dependencies.

The float32 RGBA texture contains Jr,Jtheta,Jz in A/m² and support. Manual
bilinear interpolation works without optional floating-texture filtering.
All four corners must have support, conservatively omitting boundary cells.
The ray marcher samples J(z,sqrt(x²+y²)) inside the physical texture bounds.
It never fills empty exterior, invents azimuthal variation, or evolves a state.

Front-to-back emission/absorption is an illustrative transfer function of
current (magnitude to the power 1.6), not radiation or a camera prediction. Step-length-corrected opacity
keeps appearance consistent with ray length. Threshold and opacity only affect
visibility; numeric extrema and saved arrays are untouched. Signed components
use blue/white/red; magnitude uses the displayed spectrum. Zero current is dark.

GPU errors throw to the viewer, which shows a message and falls back to the
original current samples. Texture upload happens only when the saved frame
changes. Canvas backing dimensions are capped by the shared viewer for useful
interactive performance. The original Canvas2D view remains available.

# Detailed scientific viewer shell (volume.html)

Black, stacked current-volume and magnetic-line scenes follow the supplied
scientific-render reference. Common camera controls, SI legends, actual saved
timestamps and solver-step disclosure remain visible. Current volume appears
first; the field view retains individual continuous paths. No decorative plasma
textures are used. Tube/coil geometry is an optional overlay at its true scale.

The volume builder injects DATA, VOLUME and VIEWER markers. WebGL renders into
`volume` behind the transparent `current` interaction/annotation canvas.
All controls belong to the shared viewer; the renderer owns only GPU drawing.
Opacity and visibility threshold explicitly affect display, not simulation.
Expanded source details disclose axisymmetry, interpolation, coarse physical
mesh and the difference between an illustrative transfer function and light.

This separate shell preserves all historical generated viewers. No network
libraries, fonts, accounts or external graphics are required.

Default illustrative opacity is 5; it can be reduced to inspect weak interior
current or raised to emphasize the exterior layers. No lighting is simulated.
