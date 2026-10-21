# hybrid3d_display.py

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

Volumes store little-endian RGBA16F `(Jr/s, Jtheta/s, Jz/s, support)`, where
`volume3d.value_scale = s` is the frame's peak |J| in A/m² (one for zero J).
The shader rescales its physical maximum, preventing half-float overflow.
Extrema remain in physical units; normalized half precision has about 0.05%
relative rounding. Nonfinite input is rejected.

`hybrid3d_geometry` computes flux, radius, length and density modes from the
exact displayed snapshot. Loss counters are added only for coincident history
times. The gold contour describes the azimuthal average, not 3D topology.
Ions have no persistent identity and must not be connected into trails. Their
velocity-spread equivalent includes bulk motion and is not an interior Ti.
Ion mass is taken from setup or the source gas, including D2.
The ion-dot legend discloses the speed color scale and saturation. Losses count
the entire loaded gas, including ions outside the initial separatrix.

The browser switches model labels, component availability and rendering notes
with each run. Display interpolation adds no simulation resolution or stability
evidence. Simulation/output files are kept outside version control.
