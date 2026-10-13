# Trajectory display renderer

Offline WebGL2 view of the independently reconstructed saved vectors. Derived
from the earlier showcase in a separate directory; historical assets are intact.
The current volume revolves |J| about the cylinder axis; full 3D magnetic paths
already contain the actual azimuthal field. No invented current rings.

Physical aspect ratio is preserved. The cylinder wire outline is a geometric
guide, not a conductor. Ray samples are clipped to the cylinder. The thin
half-pixel strip at its edge clamps to the last display sample; it is not vacuum.
The missing quarter-sector is only a display cutaway. Linear texture sampling
and optical opacity are presentation choices, not radiation physics.

`FieldView.setFrame` loads one actual saved state. `draw` shares camera settings
across panels and accepts export sizes. `state.localScale` selects shared or
per-frame current normalization explicitly; the controller updates its legend.
The B color range remains shared over all times. Magnetic direction and scale
are never exaggerated. Zero-current pixels have zero volume opacity.

Three-dimensional presentation of an axisymmetric calculation does not imply
three-dimensional simulation or stability. All scripts run offline.
