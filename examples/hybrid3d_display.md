# hybrid3d_display.py

When `b_response` is saved, current is its curl/μ0. Field lines still use
total B. This keeps the known external coil current out of plasma current
volumes, including when coarse Cartesian stencils overlap a nearby coil.
New runs also identify matching current/E edge support and gas-volume node
normalization; displayed current uses that exact support and normalization.
Historical prototypes retain their recorded interpretation with an explicit
notice that their boundary-current treatment has since been corrected.
The status states the available saved-time interval and, for partial exports,
the requested target and in-progress state. The timeline uses only actual
frames; no padding or relabeling makes a short run appear to cover 40 µs.
Metal-backed runs receive an explicit mixed-precision development notice,
because matched CPU/GPU shots show local-field sensitivity despite close
global budgets. Exporting saved fields does not validate their accuracy.

Builds the shared formation viewer from actual 3D hybrid snapshots:
`python -m examples.hybrid3d_display OUT RUN[::Label] ...`.
The run selector retains other existing 2D and 3D entries. Only supplied run
entries are replaced; the 2D overview is not recomputed.

- `trace` integrates 3D field lines with midpoint steps and trilinear B.
  They are not rotated copies or particle trajectories.
- `run_display_3d` exports current volumes, sparse current fallback samples,
  independently sampled ion dots and n=0 flux contours. Current is curl B/mu0,
  transformed to cylindrical components at supported nodes.
- `main` updates manifests, `runs.json` and the common viewer shell.
  Concurrent hybrid exports take a local file lock, preserving every run
  entry and a consistent viewer when multiple long simulations finish.
- `refresh_viewer` rebuilds only that shell, retaining all scientific frames.

Annular runs retain the inner tube, all coil groups, annular n=0 flux levels
and vacuum topology in the source data. Their description distinguishes
prescribed saved bank-current waveforms, coupled banks and fixed conductor
flux. Only the gold n=0
diagnostic is rotationally copied; fields and current use full 3D samples.
The particle slider changes only the displayed subset of saved random dots.
New snapshots carry actual electron pressure `pe`; their mean Te is computed
with node gas fractions and density, including thermal/heating evolution.
Conservative energy runs identify their pressure-work and boundary-loss
closure, first-order upwind diffusion, and remaining grid/loss-model limits.
They do not inherit the unresolved-energy notice of invalidated entropy runs.
Older files retain their recorded isothermal value. Coupled-bank, prescribed
waveform and fixed-flux drives receive distinct descriptions. The gold
checkbox is labeled n=0 flux outline for hybrid data.

Volumes store little-endian RGBA16F `(Jr/s, Jtheta/s, Jz/s, support)`, where
`volume3d.value_scale = s` is the frame's peak |J| in A/m² (one for zero J).
The shader rescales its physical maximum, preventing half-float overflow.
Extrema remain in physical units; normalized half precision has about 0.05%
relative rounding. Nonfinite input is rejected.
Runs with `INVALIDATED.json` are rejected before export.

`hybrid3d_geometry` computes flux, radius, length and density modes from the
exact displayed snapshot. Loss counters are added only for coincident history
times. The gold contour describes the azimuthal average, not 3D topology.
A closed pocket may be off the geometric midplane, with no midplane radius.
`flux_summary` reports that explicitly rather than inventing a radius or
requiring a centred configuration. Full saved-frame exports cover such states.
The displayed maximum radius comes from the same selected n=0 region's
`max_separatrix_m`. The full axial profile may contain infinite sentinels for
separate wall-reaching/open flux; these are not a pocket radius. Scientific
frame/manifest JSON rejects nonfinite numbers so a malformed diagnostic
cannot silently leave an unreadable browser export.
Ions have no persistent identity and must not be connected into trails. Their
velocity-spread equivalent includes bulk motion and is not an interior Ti.
Ion mass is taken from setup or the source gas, including D2.
The ion-dot legend discloses the speed color scale and saturation. Losses count
the entire loaded gas, including ions outside the initial separatrix.

The browser switches model labels, component availability and rendering notes
with each run. Display interpolation adds no simulation resolution or stability
evidence. Simulation/output files are kept outside version control.
