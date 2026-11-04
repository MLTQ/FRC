# iec_checkpoint.py

Complete full-step restart state for the existing 3D electrostatic integrator.
This changes persistence only, not forces, sources, controller or timesteps.

`write_checkpoint` atomically replaces one compressed checkpoint after the
second kick and saved diagnostics. It stores full particles, potential/charge,
electric field, electrode charges, all accumulated energy/charge/loss counters,
initial energy, exact time/step, NumPy generator state, controller internals,
history length and compute-time offset. Config, mesh digest and physics-module
hashes, numerical-library versions and architecture bind the restart. Saved
history is included, so persistence does not depend on another history file.
A failed temporary write leaves the previous checkpoint.

`load_checkpoint` rejects changed config, geometry, physical source or nonfinite
arrays. Frame/display samples are never restart state. Existing histories must
match the checkpoint time/row count; no silent truncation or overwrite follows.
The initial legacy archives need a separate measured reconstruction/replay before
they become complete checkpoints. This module does not claim they already are.

The module hashes cover the field, particle, mesh, vacuum, drive and diagnostic
modules; they do not independently guard the integrator driver. Preserve the
archived driver/source provenance as well when continuing. The Oct 25 recovery
uses a frozen driver and an exact saved-interval replay; changed drivers need
equivalent continuation evidence before adoption.
