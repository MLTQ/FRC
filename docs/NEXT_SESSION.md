# Resume here next week

Research paused at the user's request on 2026-09-30 to conserve credits.
Do not launch further studies until the user resumes. The overarching FRC goal
is open, not complete. This is a local repository with no configured Git remote;
issue state and research checkpoints are saved locally.

## What is verified

Commit `267ae93` is the last fully verified solver checkpoint: physical
regular-fluid/current/magnetic coupling, 402 full-suite tests and 16 audited
coupling cases. Frozen calculation source:
`8c9339a6eaf0c88b88bffe108d5ca2c535a8bfbf33fa4bd11b0c0aaebba6c784`.
See `REGULAR_COUPLING_NUMERICS.md`, `REGULAR_COUPLING_VERIFICATION.md` and
`outputs/regular_coupling_verification/`. This verifies an instantaneous
component, not complete moving-species evolution or FRC formation.

## Ordered work in beads

1. **frc-7mm.3 — finish and audit the current-transport candidate.** The new
   `plasma_bench/physical_derivatives.py` and `current_transport.py`, independent
   `examples/current_transport_reference.py`, study and tests are saved as an
   unfinished candidate. Four focused tests passed in 1.66 seconds; the full
   suite has NOT been rerun with these changes. The 16-case transport study is
   NOT frozen or audited. `run_case()` works, but the study CLI imports the
   **not-yet-written `examples/current_transport_records.py`**. First read the
   companion docs and `REGULAR_CURRENT_TRANSPORT_DESIGN.md`, implement the
   records/provenance and independent audit, then freeze source and run the
   degree/resolution/quadrature/density protocol including 32 elements. Check
   local, first-element and maximum errors, true continuity and current jumps,
   omitted-face negative controls, energy, axial force and angular torque.
   Preserve the initial failed arbitrary 5x negative-control assertion and its
   correction in `outputs/current_transport_diagnostics/`; no solver was changed
   to obtain that focused-test pass. Full regression and independent physical
   integration are required before adoption.
2. **frc-7mm.4 — complete species pressure/current-entropy coupling.** Work with
   the actual regular n,q,l,mz,Ce,Ci state, separate species temperatures and
   enthalpy/Poynting fluxes. Verify original parent-equation rates and unforced
   coupled evolution; do not use energy resets or fitted cancellation. Continue
   the existing fluid order/spectral limitations under **frc-7mm.1** and coupling
   work under **frc-7mm.2**.
3. **frc-7mm.5 — full moving-plasma and formation validation.** After the above,
   validate trajectories, heated moving exact vacuum, physical transport and
   ionization, ends/coils and converged formation. Preserve exactly zero plasma
   outside; no conducting/massive atmosphere. Hardware is prescribed waveforms.
   Closed flux or an axisymmetric equilibrium does not establish 3-D stability
   or particle confinement.

## Existing remote study — inspect before doing anything

The older frozen 3f study under **frc-5b5** was last observed at
2026-09-30T18:06:01.201528Z: 13/14 cases complete and audited, final nr32 contact
2050/4116 accepted steps. Host `m@192.168.0.202`, directory
`/home/m/Code/frc-tool-comparison/face-radial-20260930`, parent 147751 and workers
147839/147845/147852. This is CPU NumPy/SciPy on the GPU machine; no GPU speedup
is claimed. Status is historical, not a fresh liveness check.

No stop or restart was issued for the pause. On resumption, inspect terminal
records, process command/cwd and saved checkpoints first. Never restart a live
job or infer failure from an observation timeout. Retrieve and audit the final
case if finished. The user has already authorized transfers to this machine.
Saved observation: `outputs/regular_coupling_diagnostics/remote_status.json`.

## Updated presentation

`outputs/frc_showcase_checkpoint/index.html` is a self-contained 3-D viewer;
`poster.png` and `orbit.mp4` are the still and camera-orbit exports. Rebuild/serve
with `.venv/bin/python -m examples.showcase_checkpoint --serve --port 8779`, then
open http://127.0.0.1:8779/. This presentation uses frozen verified coupling
case 07, not the unfinished transport candidate. It samples all three current
and magnetic components, includes toroidal magnetic field in the streamlines,
and preserves the 35 mm radius / 240 mm length proportions. Display interpolation,
opacity and streamline seed spacing are illustrative. The animation changes
only the camera, not the plasma state. The older initialized-FRC showcase and
the user's existing poster remain available unchanged by this update.
