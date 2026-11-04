# hybrid3d_frc.py

Driven initialization maps the plasma-response flux separately from the
known vacuum coil basis. A t=0 start has exactly zero response field;
warm starts retain their mapped plasma current without coil-stencil current.

Driven snapshots retain both total B and `b_response` (total minus the known
coil field), plus `source_current_a`; unit coil fields are stored once in
`coil_basis.npz`. Hall force, Ohmic heat and displayed plasma current use
the response curl, so an unresolved coil stencil is not counted as plasma.

`thermal=1` activates conservative electron internal energy with Spitzer
resistivity and Ohmic heating (default gamma=5/3), replacing the legacy
fixed/bare barotropic closure and invalidated entropy transport. Snapshots retain actual node electron pressure `pe`; new
snapshots also store it for legacy closures. Diagnostics count Ohmic heat,
cutoff-removed heat and maximum electron transport fraction. See
`plasma_bench/hybrid3d_energy_electrons.md` for numerical diffusion and missing
physics. Diagnostics distinguish outgoing enthalpy, energy left on emptied
nodes, the grid energy-equation residual, explicit pressure time-centering
error and roundoff. Coupled runs also measure response magnetic energy and
the actual load-energy residual including every declared loss. This budget
is not enforced by rescaling any state.
`load_energy_work_j` is the directly integrated coil-source work, not the
bank's image work alone. A later-connected section has a nonzero flux
reference; the circuit snapshots record that exchange and waveform lag.

Starts a formed 2D FRC in the kinetic-ion/fluid-electron 3D device model.
`python -m examples.hybrid3d_frc OUT RUN t_us=10 t_end_us=6` runs a 6 µs
continuation. `dx`, `dz`, `ppc`, `floor`, `seed`, `every_us` and `snap_us`
control resolution, sampling and the declared Ohm-law density floor.
`ppc` must be a positive perfect cube (1, 8, 27, 64, ...), because positions
use an equally spaced tensor-product lattice in each gas cell. Invalid counts
are rejected instead of silently rounding: historical `ppc=16` and `ppc=32`
both loaded 27 samples/cell. Such archives must be interpreted using their
actual particle counts, not the nominal setting. New `particle_loading` setup
metadata gives the lattice shape, samples per gas cell, gas-cell count and
number of zero-weight samples removed. `particles` counts retained ions;
an initially zero-density sample does not become a numerical atmosphere.
The accepted cubic loads preserve their original positions, weights and
random velocity sequence. This validation does not alter frozen campaigns.
Actual coupled-bank 0–0.1 µs replays at 8 and 27 samples/cell retain identical
saved total B, n, pressure and ion samples versus the prior runs. The derived
response-field subtraction differs by at most 1.1e−16 T. A requested 32 is
rejected before a physical setup is written. Evidence:
`outputs/hybrid3d_loading_contract/comparison.json`.
For an actual full shot, use `t_us=0 t_end_us=40 drive=2 thermal=1`.
The duration is relative to the initial time; changing only the duration
does not turn a warm start into an initial-fill experiment. The completed
six-case v4 campaign uses 20×20×40 mm, eight ions/cell, 1 ns steps,
5% density cutoff and 0.1 µs snapshots. Its numerical heat removal and early
resolution sensitivity preclude validated configuration/confinement claims;
see `docs/HYBRID_3D.md` and `outputs/hybrid3d_full_shots_v4/summary.json`.

`gpu=1 thermal=1` explicitly selects the Metal transport prototype; default
`gpu=0` remains CPU64. `setup.compute_backend` declares the mixed-precision
path. GPU calls include explicit synchronization/casts in their wall-time
counter; full-shot timing must still include all remaining CPU work. This
does not switch the source-frozen running campaign. See
`plasma_bench/hybrid3d_metal_energy.md` for precision and conservation contracts.
`gpu=2` also enables normalized float32 atomic particle deposition, while
particle integration remains CPU64; its gather/deposit adjoint error must be
measured separately. Both modes are experimental and leave the default alone.
`arrays=1 thermal=1` optionally fuses CPU64 energy stencils and reuses their
scratch arrays, on either compute backend. The default remains `arrays=0`.
`build(..., fused_arrays=True)` is the equivalent Python API. The setup's
thermal options identify the choice. Equation, loss model and vacuum are
unchanged; existing frozen shots do not load these edits.

- `full_flux` reconstructs the source run's full flux and loads its fluid state.
  `t_us=0` reconstructs the declared initial fill/PI profile, zero bulk flow
  and current-free bias field. Older archives begin at 0.1 µs; that evolved
  snapshot is never relabeled as the initial state.
  Annular sources retain the absorbing inner wall and exact inner vacuum.
  Sources containing toroidal field/rotation are still rejected rather than
  silently dropping those fields.
- `build` maps its vector potential to Yee edges, takes the discrete curl,
  and loads thermal ions. Electrons are isothermal unless the thermal or
  barotropic closure is selected; the conductor holds flux
  and the device kernel projects the exterior gap to current-free vacuum.
  An annular vacuum also retains its meridional circulation harmonic.
- `drive=1` retains auxiliary-coil currents and follows saved vacuum-equivalent
  bank currents, using `hybrid3d_sources`. These currents are prescribed;
  3D plasma feedback on the electrical circuit is not yet included. Only
  the coil waveforms are replayed, never the saved 2D plasma response. Runs
  cannot extend beyond those waveforms. `drive.json` stores the actual input.
  Auxiliary-coil sources require this option.
- `drive=2` starts the full switched-bank circuit and plasma together at t=0,
  with 3D current reciprocity supplying the bank load (`hybrid3d_circuit`).
  This new coupling needs independent linkage/convergence verification.
- `gamma_e>1` selects the existing barotropic adiabatic electron closure;
  otherwise electrons remain isothermal. Neither closure includes Ohmic
  electron heating or ionization. `dt_ns` can impose a smaller fixed ion step
  for formation, where B increases substantially after initialization.
- `diagnostics` combines `hybrid3d_geometry` n=0 flux/modes with particle,
  loss and energy counters. `twod` rows are nearest saved reference times,
  explicitly retained for comparison.
  Particle losses cover all loaded gas, not a cohort selected inside the FRC.
  Diagnostics include mean electron temperature, adiabatic internal energy,
  divergence, ion Ωdt and cell-travel bounds. Circuit residual is reported
  separately and is not a closed plasma+boundary energy budget.
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
Thermal setup also records the0.8 electron transport loss CFL. Histories count
energy substeps and the maximum effective outgoing-loss fraction, including
advection, wall enthalpy and pressure expansion. Fields retain their own RK
bound; this numerical correction does not rescale U or change the loss model.
They also identify the matched plasma-current/E edge support and gas-volume
node normalization. Earlier full-shot trajectories with unmatched conductor
sheet currents are invalidated; their artificial Hall power is not a mode.
Snapshots are atomically renamed into place, so a concurrent viewer cannot
read a half-written archive. A terminal snapshot and diagnostic are retained
even when the duration is not a multiple of the sampling interval.
