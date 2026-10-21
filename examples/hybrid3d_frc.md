# hybrid3d_frc.py

Starts a formed 2D FRC in the kinetic-ion/fluid-electron 3D device model.
`python -m examples.hybrid3d_frc OUT RUN t_us=10 t_end_us=6` runs a 6 µs
continuation. `dx`, `dz`, `ppc`, `floor`, `seed`, `every_us` and `snap_us`
control resolution, sampling and the declared Ohm-law density floor.

- `full_flux` reconstructs the source run's full flux and loads its fluid state.
  Annular sources and sources containing toroidal field/rotation are rejected:
  their inner wall and extra fields have not yet been mapped into this driver.
- `build` maps its vector potential to Yee edges, takes the discrete curl,
  and loads thermal ions. Electrons are isothermal; the conductor holds flux
  and the device kernel projects the exterior gap to current-free vacuum.
- `diagnostics` combines `hybrid3d_geometry` n=0 flux/modes with particle,
  loss and energy counters. `twod` rows are nearest saved reference times,
  explicitly retained for comparison.
  Particle losses cover all loaded gas, not a cohort selected inside the FRC.
- `save_snapshot` stores float32 face B, node n, and an independent random
  ion-position/velocity sample. Sample identities are not persistent: these
  snapshots cannot be joined into particle trails.
- `save_history` atomically replaces the scalar history, including its initial
  and terminal diagnostics.

`setup.json` describes the mapping; `history.json` holds scalar diagnostics;
`snap_NNNN.npz` feeds `hybrid3d_display`. Reported n=0 closed flux does not
establish 3D stability. The electron closure and absorbing 3D boundaries differ
from the 2D formation model, so their continuations need not coincide.

New runs record ion mass, seed, requested duration, sampling intervals,
actual cell spacing and the complete `HybridParams` electron/Ohm-law settings.
Snapshots are atomically renamed into place, so a concurrent viewer cannot
read a half-written archive. A terminal snapshot and diagnostic are retained
even when the duration is not a multiple of the sampling interval.
