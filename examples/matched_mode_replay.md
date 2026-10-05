# Every accepted state and stage replay

Reruns a complete successful signed trajectory with exactly its original mesh,
signed amplitude, duration, step count and tolerance, selecting the explicitly
declared guarded or unguarded solver. Existing directories cannot be overwritten.
The separate generic archive verifier requires identical frozen source/runtime,
every raw file/key/value (including every accepted stage and endpoint), and
all non-runtime summary/iteration metrics. Physical reconstruction fingerprints
are verified separately by repeating the independent audit.

A failed rerun records its last accepted state and cannot produce successful
verification. Replay demonstrates reproducibility, not physical convergence.
