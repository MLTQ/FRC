# Smooth coupled provenance

Hashes the complete local solver/oracle/measurement import closure, package
initializer and provenance code. Independent of the historical 21/18 baseline
identity. Each measurement freezes the complete identity at start/end and saves
full state arrays, inputs and NPZ hash. Report-only helpers/tests need a separate
snapshot. A source edit during a live run invalidates that run rather than
silently relabeling it under a newer method.
