# Immutable Radau protocol and audit

Freezes explicit source closure/runtime before each single trajectory and
refuses overwrite. Accepted-state checkpoints and terminal raw arrays use
hashes; exceptions retain the actual last accepted checkpoint and traceback.
Only a completed trajectory receives case.json. Sources must stay unchanged
throughout the calculation.

The audit action verifies original source/runtime/raw hashes, independently
reconstructs the contact/wave fields and continuous physical energy at both
quadrature orders 16/20, and optionally repeats the complete trajectory. Every
saved array and non-runtime metric must replay bitwise. No terminal-success
summary is written on a failed check. Reuses the independently documented
implicit_plasma_audit.measure formulas; this is separate from the production
grid work. Record numerical temporal damping and convergence explicitly.
