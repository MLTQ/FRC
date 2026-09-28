# Active objective: a credible experimental FRC simulation

The goal is a numerically verified and physically justified simulation in which
we can vary geometry, plasma and prescribed coil waveforms to investigate FRC
formation. Passing the existing stationary-conductor tests is not completion.

## Completion evidence still required

1. Conservative axisymmetric compressible dynamics: mass, momentum, pressure/
   total energy, magnetic advection and resistivity, with compatible div B.
2. Verified plasma–vacuum coupling and suitable chamber/end/outer boundaries;
   no undisclosed conducting or massive numerical atmosphere.
3. A defensible transport/ionization model for the selected regime, with local
   scale/collisionality checks; Hall/two-fluid or kinetic treatment if required.
4. Independent tests of waves, shocks, magnetic forces, diffusion/heating and
   equilibrium preservation, including convergence and energy budgets.
5. Converged formation experiments with separate reversal, topology, force
   balance, lifetime and model-validity criteria. Robustness to time step, mesh,
   boundary and transport choices; no claim of stability from axisymmetry alone.
6. Scriptable and interactive experiments with complete reproducibility records
   and diagnostics that do not overstate the supported conclusions.

## Current evidence

The initial turn was progress: Stage 0/1 were implemented with 33 passing tests
and three waveforms plus a 27-case sweep. Their topology durations are not
converged and their illustrative argon parameters raise missing-physics warnings.

The current goal turn is progress: a conservative dynamic r-z MHD kernel now
passes independent cylindrical fast-wave, 2-D flux-advection, exact Sod shock,
resistive-decay/heating, axis-balance, divergence and temporal-order checks.
The full suite now has 44 passing tests. Direct analytic J×B acceleration and
a static magnetic-pressure equilibrium expose the remaining first-order axis
error and artificial equilibrium flow; both are explicitly measured.
Measured convergence and conservation are in `docs/MHD_VERIFICATION.md` and
`outputs/mhd_verification/metrics.json`. It is a bounded verification solver,
not the coupled coil-driven formation experiment.

## Current development step

The first-order Rusanov/UCT reference kernel is in place. Its spatial errors
decrease under refinement, but its numerical diffusion is too large to ignore
in topology-sensitive formation studies. Next, improve the spatial method and
verify magnetic equilibrium/force balance, then connect it to the independent
electromagnetic/vacuum subsystem with a measured source-work budget. Do not
simply enable `Experiment(stage=2)` with periodic/metal-wall boundaries.

One physically relevant coupling route is a pre-ionized gas filling the actual
chamber, with fluid only inside the chamber and a true magnetic vacuum outside.
That avoids adding an artificial massive exterior. It still requires validated
quartz/end boundary treatment and must address rarefaction, changing conductivity
and formation of an internal plasma–vacuum interface if a detached torus forms.
It is an implementation candidate, not a completed or validated model.

## Completion audit at this checkpoint

| Requirement | Authoritative evidence | Status |
|---|---|---|
| Conservative dynamic core | mhd modules, independent tests, MHD verification report | Implemented for stated bounded domain; first-order spatial accuracy |
| Physical source/plasma/vacuum coupling | Original Stage 2 guard still raises; beads frc-vacuum | Incomplete |
| Regime-appropriate transport / Hall physics | Fixed eta/sigma only; beads frc-transport | Incomplete |
| Full independent numerical verification | Wave, shock, diffusion/heating, advection, direct J×B and straight-field equilibrium data; no FRC equilibrium/strong 2-D formation validation yet | Partial |
| Converged FRC formation | Only older stationary-conductor results, unconverged topology duration | Not demonstrated |
| Reproducible experimental interface | Stage 0/1 CLI/playback; dynamic benchmark scripts and metrics | Partial; no coupled dynamic experiment UI |

Beads issues `frc-goal`, `frc-mhd`, `frc-vacuum`, `frc-transport`, and
`frc-formation` preserve the remaining scope. No external process is being
waited on. The optional request for target plasma parameters is not a blocker
to numerical verification or development.

The full objective remains active until all completion evidence is present.
