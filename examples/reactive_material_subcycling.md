# Full charged argon with explicit reaction intervals

Retains the full20ns neutral/ion/electron ground-ladder argon problem and its
physical masses,1e14/m³ initial component densities,10mm half-width and
0.1/0.2/2eV temperatures. The precise geometry and vector-assembled source are
unchanged. Physical-width coarsening is absent. `--coalesce-ulps` now defaults
to zero, retaining every distinct density cut. The earlier recorded10ns/5ns
attempts used8-ULP internal cleanup; that remains an explicit comparison option
and still protects every species endpoint. Existing atomic/closure limits remain.

`--source-step-ns` chooses a disclosed source interval (default10ns), not the
motion timestep. Within it the material equations/fields use actual RK4 steps
at half their current acoustic/plasma bound. Every accepted internal dt and
physical budget is retained separately. Reaction states can change rapidly
through the motion, so adequate source timing must be demonstrated by complete
10/5/2.5ns or finer comparisons; no convergence claim follows from completing
the first chosen interval. No implicit transport or enlarged physical mass.

Each completed coupled interval saves the full precise state and reservoir
budgets. Internal transport stages save times/budgets, not full geometry
snapshots; they are explicitly split-stage records, not completed coupled states.
Source/input identities, initial state and actual internal step sizes retain
reproducible evolution. A rejected source/transport interval retains its last
source and internal diagnostics but does not advance completed coupled time.
No clipping, hidden retry or shortened success. Disk reserve remains256MiB.

Use `python -m examples.reactive_material_subcycling FRESH_OUTPUT
--source-step-ns 10 --cells 32`. This is an opt-in coupled research calculation,
not a replacement for the required spatial, source-time and transport-time
convergence, physical collision/field/wall integration or FRC formation evidence.
