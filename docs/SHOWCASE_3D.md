# Current density and magnetic flux density in 3D

The user-requested showcase is a presentation of the existing completed
`source_g1_32` moving resistive verification run. It is separate from the latest
prescribed-flow bulk thermal study. No new physical simulation is inferred from
the rendering. The run starts with reversed flux already initialized; it has
unresolved heated-edge and spatial accuracy documented in SURFACE_VERIFICATION.

- Interactive offline artifact: `outputs/frc_showcase/index.html`.
- Exported and visually inspected 2800 × 1900 PNG: `outputs/frc_showcase/poster.png`.
- Data, source identity, physical parameters and sampling audit:
  `outputs/frc_showcase/data_audit.json`.
- Snapshot times are the actual saved 0, 20 and 40 ms states. Colors share fixed
  limits of 0–11.80436 kA/m² and 0–1.659691 mT across all times.

The current-density volume and toroidal current rings use the archived weak
current. The sampler inverts each actual deformed quadrilateral and evaluates
its original Q1 field, rather than smoothing the simulation into an analytic
torus. Current quadrature replay errors at the three times are 1.24e-16,
7.01e-16 and 4.15e-16 relative L2. Current is zero outside the material in the
source and sampled data. Display texture interpolation can soften a boundary
over a pixel; it is not an added atmosphere.

Magnetic field is evaluated from the archived flux coefficients and the actual
deformed mesh gradient, including the magnetic vacuum. Display field lines are
contours of the sampled poloidal flux, copied into meridional planes around the
axis. Their maximum flux discrepancy relative to peak flux is 8.66e-4,
1.46e-3 and 1.29e-3. Tube thickness, line spacing, optical opacity and glow are
illustrative. There is no invented helical field or 3D instability calculation.

The viewer has synchronized drag rotation, bounded zoom, camera auto-orbit,
three actual saved time selections, a display-only cutaway and a labeled poster
export. Browser inspection verified all three snapshots with fixed legend
limits, cutaway, orbit, synchronized drag, reset and successful PNG export.
The poster was opened and visually checked for readable labels and unclipped
graphs. The standalone HTML embeds all data and scripts; no CDN or external
fonts are required.

The localhost preview is currently served on port 8778 by session 11395 and its
browser tab is retained as a user deliverable. It is a preview process, not an
active scientific calculation. The saved HTML also works after that server is
stopped; its export then uses a browser download.

Rebuild data with `python -m examples.showcase3d_data`; assemble with
`python -m examples.showcase3d`. Optional local preview:
`python -m examples.showcase3d --serve --port 8778`.
