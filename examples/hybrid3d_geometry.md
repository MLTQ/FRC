# hybrid3d_geometry.py

Read-only geometry shared by the hybrid run driver and its saved-state viewer.
It preserves the driver's original annular averaging and midplane definitions,
so viewer diagnostics come from the displayed snapshot rather than a nearby
history row.

- `poloidal_average`: integrates annular Bz face averages into psi and calls
  `closed_region` inside a solid tube or `annular_closed_region` outside an
  inner tube. The latter excludes flux touching either radial wall. Returns
  r, z, psi and the n=0 region; radial boundaries use the next outer bins.
- `density_modes`: returns instantaneous midplane n=1,2,3 amplitudes in the
  annulus from 0.3 to 0.9 wall radii. Particle noise is included.
  An inner tube is excluded from this diagnostic annulus.

Both consumers supply the saved Cartesian box/device geometry and correctly
staggered Bz or node density. These quantities describe the axisymmetric
component; they do not prove full 3D topology or instability growth.
