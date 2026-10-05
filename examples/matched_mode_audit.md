# Signed full-trajectory continuum and conservation audit

Requires complete successful opposite-amplitude trajectories with identical
meshes, duration, steps, tolerances, reference and frozen source. Verifies every
saved-step and initial/final raw hash. Independent root-polynomial quadrature
16/20 compares centered physical fields, both pressure responses, full E and J
with the separately derived continuum eigenmode times exp(lambda t). Retains
initial representation, axis/wall/max error and component-wise projected growth.
Projection measures amplitude; it never changes or fits a simulation state.

Each sign has a separate original continuous energy/particle/entropy audit,
including magnetic bias cross work and positivity. Small-amplitude and time/mesh
controls are still required; good energy or growth alone is not physical field
accuracy. Full solved electric fields and strict induction diagnostics are
retained. Expanded physical arrays have dtype/shape/value hashes and can be
reconstructed exactly from the complete saved coefficient arrays.

Freezes a separate complete source closure and writes terminal only after all
checks pass. This audit does not itself replay evolution. Full bitwise replay
of every stage/state and metric is a separate gate. The fixture remains bounded
periodic fully ionized plasma, without formation, vacuum motion or end losses.

The fixture returns seven entries; only the reference, fluid, system and bias
are needed for reconstruction. Initial unchanged-state controls exercise the
complete audit path before any trajectory result is frozen.
