# Original versus streamed complete-workflow comparison

Requires two completed successful signed physical workflows with identical
mesh, reference hashes, amplitudes, actual solver diagnostics and all raw
states/rates/face-integral arrays. Rechecks input hashes and frozen source
closures. Every coefficient must agree bitwise; no tolerance is applied to the
solver or its physical diagnostics.

Reports every quadrature norm's absolute and relative change, plus physical
field fingerprint equality. Diagnostic floating-point reduction/batching
differences are measured, not hidden with an equality assertion or a loosened
physical acceptance limit. Near-zero diagnostics retain both scales. A large
norm discrepancy requires investigation even if the full states agree.

run freezes inputs and source, refuses overwrite and can itself be replayed
with wall_diagnostic_replay. It also accepts two batching sizes for a complete
workflow control. This is end-to-end evidence, not a unit test or accuracy
certification for the underlying physical reference.
