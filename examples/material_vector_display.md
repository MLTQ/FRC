# Actual moving full-vector 3-D display

Reads a completed successful trajectory and reconstructs five actual canonical
states with its declared field order, quadrature, resistance and species model.
Current is en*w from the physical compatible Ampere solve. No time interpolation,
rescaling of scientific states or further simulation is performed.

`trace` follows the full interpolated B vector in Cartesian arclength, so it
does not require positive Bz or force an open field line. It stops outside the
sample support, near weak field, near the unresolved axis, after a displayed
loop or at the fixed display length. Linear physical-space interpolation and
finite tracing steps are presentation approximations, not topology certificates.

The current cloud samples actual quadrature values and revolves them through
16 azimuths. Magnitude and all signed cylindrical components are retained. The
outline comes from actual moving boundary nodes. One physical length scale is
used in all directions, with limits covering the whole moving outline. A PNG/PDF
shows the last state's field lines and current magnitude; the offline HTML exposes
actual-time selection, linked rotation/zoom and explicit color-scale controls.

The source record retains raw-state and summary hashes without replaying unrelated
scientific campaigns. The visible scope reports coarse initial B/J errors,
initialized flux, absent coils, upwind diffusion and unresolved physical validity.
No field lines are claimed to be particle paths or 3-D stability evidence.
