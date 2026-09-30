# Archived-field 3D presentation sampling

Reads the completed `source_g1_32` resistive moving fixture and verifies its raw
archive hash. No simulation is rerun or changed. The deformed material and
magnetic-vacuum Q1 coordinates are recovered from the recorded configuration.
Triangles locate the containing quadrilateral; Newton inversion then evaluates
the original bilinear field, its physical magnetic gradient and weak current.
Current is exactly zero outside material support. Saved current quadrature is
replayed independently of the display grid before export.

The display domain stays inside the ordinary magnetic mesh wherever values are
sampled; points outside its support are zero and masked before contouring.
Entirely unsupported sample batches return an empty support without reduction.
Contours of the sampled psi are approximate field lines, with their discrepancy
against the original Q1 flux recorded. Toroidal current rings and repeated
poloidal contours express axisymmetry, not a 3D solution or helical field.
Color limits are common to all three actual saved times. Magnitudes use SI
current and mT magnetic density; physical lengths remain meters.

`build` returns compact presentation arrays plus source identity and limits.
The command saves data and a small audit in `outputs/frc_showcase`. Interpolation
and visual opacity are presentation choices, not changed scientific states.
The unit-scale initialized FRC retains its documented unresolved heated-edge
and mesh sensitivity. This is not coil-driven formation or stability evidence.
