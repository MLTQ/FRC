# Complete reconstruction of the tighter-bound rejected input

Validates the terminal original failure, raw input and source closure, rebuilds
the actual full nonlinear fluid/current/species forcing, and observes the first
strict split-inverse attempt using the same direct-metric setup. Base, tangent
and signed input arrays must match exactly. The original first-attempt failure
message must match; all four iterate observations and raw actions are saved.

This is a complete captured-state diagnostic workflow, not an accepted full
plasma derivative or a substitute for the original five-attempt failure replay.
The result explicitly distinguishes source/input validity, matched rejection
and atomic unchanged input. Rebuild and replay the entire workflow into a new
directory and compare all arrays/non-runtime metrics. No unit tests are used.

Existing destinations reject before work. A mismatched observation is saved
as inconsistent and then rejects; it cannot be advertised as verified evidence.
All old sources and failed inputs remain unchanged. Governing equations and
unchanged fixture assumptions are in SPLIT_RESIDUAL_PATH_DESIGN.md.
