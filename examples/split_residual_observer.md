# Unchanged split-inverse iteration with observations

Exact numeric copy of SplitInverseWallCore.solve. It adds only calls to
split_residual_paths.observe after evaluating the original residual components.
All solver matrices, original residual definitions, CG/refinement budgets,
acceptance and failure messages remain unchanged. Observation arrays are never
used to propose or accept a correction.

The exact saved-input workflow must reproduce the archived first-attempt
rejection. A fresh full reconstruction/replay must match every captured array
and non-runtime metric. This is a diagnostic, not a replacement plasma solver
or a claim that an unsuccessful full derivative was accepted.
