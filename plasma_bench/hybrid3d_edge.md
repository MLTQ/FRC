# hybrid3d_edge.py: a sharp plasma–vacuum edge (design)

A prototype for bead frc-ibl.11.9. It is a subclass of `VacuumMaskHybrid` in
new files; there are no shared solver changes.
- Analytic check: `examples/edge_compression_check.md`.
- Startup replay: `examples/vacuum_mask_replay.md`.

## The defect it targets
The floor closure in `DeviceHybrid` uses ne = max(n, n_floor) in the ideal and
Hall terms of Ohm's law. Below the floor (thin nodes) this has three effects:
- **The field slips.** It moves at (n/n_floor)·u_e, not u_e, so it is left
  behind departing ions and builds current where they were.
- **Thin ions are under-pushed.** Each feels (n/n_floor) of the J × B/n it
  should, so the magnetic piston does not sweep them up.
- **Emptied gas becomes a resistor.** With no ions, Te falls to the 1 eV
  Spitzer cutoff: η/μ0 ≈ 820 m²/s, a magnetic diffusion time of about 0.5 µs
  across a 20 mm cell. The coil field has to diffuse through it. That diffusion
  time is comparable to the implosion, the worst case for Joule loss.

Measured in the v4 startup (`examples/vacuum_mask_replay.md`):
- 68 J of numerical removal by 1 µs, mostly Joule heat in emptied nodes;
- it falls only as about N^−0.25 with ions per cell;
- deleting the current when a cell empties (the plain mask) costs 94 J.

## Physical requirement
- **Where there are no resolved carriers, the field is vacuum:** J = 0 and
  curl B = 0 (Darwin limit), set at once by the boundary data. This is also the
  project rule: exact vacuum outside the plasma, never a conducting or massive
  numerical atmosphere.
- **Inside the resolved plasma, the Hall Ohm's law holds with no floor:**
  E = −u_i × B + J × B/(ne) − ∇p_e/(ne) + ηJ.
- **At the interface,** a surface current separates the two. The plasma side
  sees the plasma's own density, not one diluted by vacuum cells.

### Resolution limit
One macro-ion per cell is n = w/V: 1.5·10¹⁹ m⁻³ in the v4 shot, 1.4 times the
floor. The field stops being shielded (electron skin depth larger than the
cell) only below about 7·10¹⁶ m⁻³. Between the two, a real plasma conducts, but
8 ions per cell cannot represent it. The model has to assign this range to
vacuum or to a sub-grid model of unresolved conducting plasma. This prototype
assigns it to vacuum, as the rule requires. That choice is declared, not hidden.

## Model
1. **Classification** (inherited from `VacuumMaskHybrid`). Plasma cells are
   gas cells with at least `min_ions` macro-ions at the start or end of the ion
   step. All other cells are exact vacuum: gap-projector potential field, no
   Ohm's-law E, J = 0. Each change of the mask is projected once, and its energy
   change is counted.
2. **Interface floor.** At a node with plasma-cell fraction f_p and gas fraction
   f_g, the ideal and Hall terms use ne = max(n, n_floor/s), with
   s = min(f_g/f_p, `max_scale` = 2).
   - The cap exists because, without it, a node whose only plasma neighbour holds
     one distant macro-ion has a near-zero CIC share. The floor there then drops
     eightfold, and so does the timestep.
   - This is the same as normalizing n, n u_i and J to the plasma part of the
     node's dual volume, then applying the floor. It is the wall's gas-fraction
     normalization applied to the moving plasma boundary.
   - In interior plasma f_p = f_g, so nothing changes.
   - With `min_ions`·w/V ≥ n_floor, the floor should bind nowhere in the plasma.
     The diagnostic `floor_bound_plasma_nodes` checks this.
3. **Energy consistency.**
   - The ideal and Hall terms share ne, so the u_i·(J × B) terms still cancel in
     E·J_e, and field–ion exchange stays exact.
   - The electron-pressure term keeps the closure's own max(n, n_floor). That is
     the field `EnergyElectrons.hall_pressure_work` uses, so pressure work stays
     consistent.
   - Spitzer η uses the closure's own density.
4. **Electron edge** (`EdgeEnergyElectrons`, on by default). The same interface
   density is used in the electron-pressure force (Ohm's law), the ion pressure
   work, the Hall pressure work and the compression velocity n u/ne. These are
   one field and its adjoint, so the exchange stays exact apart from explicit
   time-centering. Two more changes:
   - No internal-energy flux crosses a face that touches a vacuum node.
   - Energy left on vacuum nodes is removed and counted.

   A first version advected with the unfloored n u/n but kept the floored
   pressure force. That broke the adjoint pairing and was replaced.
5. **Timestep.** The base whistler/Alfvén bound assumes ne ≥ n_floor. It is
   scaled up by n_floor/min(ne) whenever an interface node goes lower.
6. **Unchanged:** particle absorption, gas-fraction moments, the electron
   closure's gas dual volumes, and all counters.

## Declared approximations
- **Unresolved plasma is vacuum.** Below `min_ions` per cell it carries no
  current. Its ions are test particles in the vacuum B, pushed by the node E of
  Ohm's law. No inductive vacuum E is computed, so they miss the E × B drift
  that a rising coil field imposes.
- **Mask transitions still project.** The energy of current left in a cell
  when it becomes vacuum is lost, and counted as `mask_projection_response_j`.
- **Interface electron energy is diluted.** Electron energy at interface nodes
  still spreads over the node's gas dual volume, including its vacuum part, so
  Te there is underestimated by f_p/f_g.
- **The interface is a staircase and noisy.** It follows PIC occupancy.

## Verification plan
1. **Adiabatic θ-pinch compression** (`examples/edge_compression_check`):
   analytic radius and core field, and zero dissipation. Run for the floor,
   mask and edge models.
2. **v4 startup replay** over 0–1 µs: energy budget against the baseline, and
   ion-count sensitivity.

## Status (2026-10-21): verified on the analytic case; no gain for the FRC startup
**Analytic compression** (`examples/edge_compression_check.md`):
- Floor, mask and edge all reach the analytic compression geometry within
  1–5%.
- The floor medium deletes 182 J of 226 J of source work. Its empty annulus is a
  resistor.
- A sharp edge cuts the work to the compression's own 52 J.
- The electron edge then cuts the edge removal from 25 to 16 J, at a 0.6 J
  time-centering residual.
- The plasma-side interface floor adds nothing measurable.

**v4 startup, 0–1 µs** (`examples/vacuum_mask_replay.md`):
- The edge model removes 104/102 J (seeds 0 and 1), against 68 J for the floor
  closure. It still projects away 92–93 J of current as cells empty.
- Raising the vacuum threshold makes it worse: 121 J at 2 ions per cell, 133 J
  at 4.
- In the startup, the formation current flows in the sparse plasma. Declaring
  that plasma vacuum deletes the current; flooring it dissipates the current.
  Neither is physical.

**Verdict:**
- Keep the production floor closure for formation runs.
- Use this model and the compression check as a reference for configurations
  that start with vacuum inside the gas.
- The electron-edge flux block is the only piece that pays off in both tests.
  It belongs in the closure's own lane (Codex), if adopted.
- A real fix for the startup needs the current-carrying sparse plasma resolved,
  not reclassified. With uniform ion counts, 64 per cell still removes 45 J.
