# Electron measurement records

Hashes include all production, analytic reference, case and measurement files
used in this isolated component study, plus Python/NumPy/SciPy versions.
The fingerprint is checked before every save and across the full collection.
Lossless NPZ arrays have individual hashes; strict JSON retains configurations,
terminal or expected-failure status, histories and error norms. Zero reference
norm returns null rather than inventing a relative denominator.

Report modules have separate provenance and must never relabel measured data
after a production change. `completed` means the declared measurement finished,
not that every evolution succeeded or the formation goal is complete.
