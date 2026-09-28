# Reproducible outputs
`Result` stores named numeric arrays and JSON metadata; compressed archives load
with pickle disabled. `provenance` includes complete typed configuration,
dependency/platform versions, code SHA-256 fingerprint, units, initial state,
boundary, energy and current-time conventions. Timestamps differ between runs;
deterministic physical arrays should not. Script callbacks require sampling
before export. Output directories are explicit, never silently global.
Stage-specific boundary and electric/current time conventions distinguish vacuum
instantaneous fields from conductor backward-interval fields.
