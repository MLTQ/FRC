# Continue from the 2026-10-04 checkpoint

The user resumed on October 4. The overarching FRC goal remains open; this is
verified component progress, not a completed formation simulator. This repository
has no configured remote. Research and issue state are saved locally.

## Verified now

- Current advection, drift stress, stretching and self-drift: 16 frozen cases,
  bitwise replay, independent order-16/20 continuous-field audit. **408 full
  tests passed in 228.47 seconds.** Source `2f682eece854d54b3b0d33026709f8d34617cc4633d8b3875a0f3e2b8e0d9540`.
  The study CLI is complete; the missing records/audit modules from the prior
  pause have been implemented. frc-7mm.3 is closed. Read CURRENT_TRANSPORT_NUMERICS
  and CURRENT_TRANSPORT_VERIFICATION for exact scopes and measured limits.
- The old remote source-3f nr32 contact finished all 4116 steps and passed its
  original-platform final-step/physical audit. All 14 face-radial trajectories
  are complete and audited; frc-5b5 is closed. Evidence is in
  `outputs/face_radial_remote_completion/`. No restart occurred.
- Four species-pressure/entropy charge-compatibility probes, both quadrature
  orders replayed bitwise. Source c6a23e1a, `outputs/species_charge_probe/`.
  This diagnostic is separate from the 408-test solver baseline. It is not an
  implemented species transport operator.
- All 13 preceding calculation source closures are unchanged. Do not mutate
  their archived source definitions or rerun into existing immutable records.

## Next action: frc-7mm.4

Read `REGULAR_SPECIES_COUPLING_DESIGN.md` and the relevant module companions.
Derive the shared current/particle flux and both species entropy updates using
the actual regular n,q,l,mz,Ce,Ci state and physical h_e,h_i energy gradients.
The identity grad(p)/n=grad(gamma*h)-h*grad(s) exposes the required charge-flux
pairing. Compatible potential gradients annihilate current work, but their local
approximation error is nonzero and must converge; a global cancellation alone
is insufficient. Do not silently change density continuity or repair heat.

First implement and verify original local pressure and entropy rates with
nonuniform density, independently varying species temperatures and nonmatching
element traces. Measure axis/wall/max errors, species entropy, charge/particle
continuity, torque and physical versus grid energy. Then combine the verified
components and run unforced complete-parent wave/contact trajectories, including
enthalpy/Poynting boundaries and convergence. The lower-order continuous-energy
and earlier fluid order/spectral limitations remain recorded under frc-7mm.1.
frc-7mm.2 remains open for full coupling; frc-7mm.5 covers later adoption/formation.

## Later formation gates

Preserve exactly zero exterior plasma in heated moving-vacuum runs. Physical
transport/ionization, ends/coils and converged formation remain required.
Hardware stays abstract prescribed waveforms. Closed flux or axisymmetric
equilibrium does not establish 3-D stability or particle confinement.

## Machine and presentation

The authorized remote host is `m@192.168.0.202`; the completed old study is in
`/home/m/Code/frc-tool-comparison/face-radial-20260930`. The original PIDs were
absent on resumption and terminal evidence was audited/retrieved. No remaining
remote run is expected from that study; revalidate actual process identity
before any future action. Existing host authorization includes data transfer.
It was CPU NumPy/SciPy work, not a GPU-accelerated solver.

The six-panel convergence report is in
`outputs/current_transport_diagnostics/convergence.png` and `.pdf`.
The earlier interactive field presentation is preserved in
`outputs/frc_showcase_checkpoint/` (HTML, poster, MP4). It still depicts the
saved coupling case 07, with camera orbit only, and its historical status label
is from the previous checkpoint. Rebuild/serve via
`.venv/bin/python -m examples.showcase_checkpoint --serve --port 8779`.
The user's modified `outputs/frc_showcase/poster.png` remains outside our commits.
