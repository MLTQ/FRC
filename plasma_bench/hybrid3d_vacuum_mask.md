# hybrid3d_vacuum_mask.py

A prototype of a **dynamic vacuum mask** for the 3D device hybrid (bead
frc-ibl.11.9). It is a subclass of `DeviceHybrid` in a new file; no shared
solver module changes. Driver and results: `examples/vacuum_mask_replay.py`.

## What it changes
In each ion step, gas cells with fewer than `min_ions` macro-ions (default 1),
counted at both the start and the end of the step, are treated as **exact
vacuum for the field**:
- Their faces join the gap. B there is the potential field from `GapProjector`
  (B = −∇χ, div B = 0, J = 0), with Neumann data from the plasma-side faces and
  conductor faces.
- Edges that touch only masked or conductor cells get no Ohm's-law E, carry no
  current and get no Joule heat. With `min_ions = 1`, every edge that still
  carries current touches a cell holding an ion, so both of its CIC nodes have
  n > 0: the empty-node Joule heat cannot occur by construction.
- Circuit linkage (`CircuitDrive.linkage`) reads `model.edge_mask`, so it uses
  the same current support.

What it keeps on the true gas geometry (`types`): particle absorption at the
tube walls, the node gas-fraction density normalization, and the electron
energy closure (`EnergyElectrons` link and wall areas). Electron energy still
lives on every gas node. Thin nodes on the plasma–vacuum interface still use the
Ohm's-law density floor and the retained-heat fraction (n/n_floor)^γ, so the
cutoff sink is reduced, not removed.

## Declared costs
- **Mask-change projection loss.** Changing the mask re-projects the present
  field once. With the boundary normal B fixed, the potential field is the
  minimum-energy field in the newly masked region (Thomson's theorem), so any
  current still flowing there is removed and its field energy is lost. This is
  counted as `mask_projection_response_j` (response field) and
  `mask_projection_total_j` (total field) and must be added to the
  load-energy residual when comparing.
- **No harmonic field in vacuum pockets.** A masked region that is a ring around
  the axis, disconnected from the outer gap, cannot hold a circulating harmonic
  field; the scalar potential erases it, and the loss appears in the same
  counter. Periodic and annular devices, whose gap already needs harmonics,
  are refused.
- **The mask follows PIC occupancy, so it is noisy.** Interface cells flip as
  single macro-ions move. The union of start and end occupancy removes the
  one-step lag but not the noise.
- The projector is rebuilt (sparse LU) whenever the mask changes; the wall time
  is counted in `mask_build_wall_s`.

## Status
A bounded experiment for the 0–1 µs startup, not a validated vacuum model. A
physical treatment would track the plasma–vacuum boundary with resolved
carriers (or a kinetic electron model). This prototype only tests whether the
numerical resistive medium, not physics, sets the low-density sinks.
