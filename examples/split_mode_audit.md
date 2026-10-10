# End-to-end physical and conservation audit

Reuses the full independently reconstructed matched-mode audit with the actual
split-inverse system for full E and equation diagnostics. Quadrature16/20,
both species pressures, all physical fields, projected growth, original energy,
particles, entropy, positivity, complete input-step hashes and terminal/source
checks remain unchanged. The benchmark is radial1,N64, not radial16,N128.

`run` audits a complete opposite-amplitude pair; incomplete or failed runs
cannot pass. `compare` compares independently audited states on identical
meshes and durations, using split-inverse full E consistently for both. A
comparison with an older audit therefore measures the actual saved state
change, not the older solver's different E reconstruction. It does not claim
bitwise equality between different solver sources or prove time order.

Separate source closure records the opt-in dispatch. Temporary serial module
bindings restore on success or failure. Run in a separate process for parallel
variants. Verify the complete audit and comparison by fresh replay, without
unit tests, in accordance with the user's preference. No prior source changes.
