# Immutable single implicit trajectory

Freezes the complete explicit calculation import closure and runtime, records
actual accepted progress with copied state arrays and their hashes, then saves
the final complete arrays/metrics or terminal exception. Existing directories
are never overwritten. Source identity must stay fixed during the run.
Progress records are observations, not proof a process remains live or a
trajectory completed. A separate audit must replay and independently integrate
the terminal accepted states. CLI defaults cover one full contact/wave period.
`--method` selects midpoint or two-stage Gauss; the last accepted pre-failure
state is retained when the nonlinear solver rejects a later step.
