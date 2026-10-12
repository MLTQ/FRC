# Guard complete scheduled-workspace refinement

Reuses the unchanged F64/256, radial16,N128, epsilon1e-4, tolerance1e-10
physical protocol with selected field nr/nz. Only the scheduled-workspace
driver changes. The control archive explicitly records that source selection.

Retains8GiB own-worker peak RSS,256MiB free disk and5s bounded own-process
termination after SIGINT. Fresh output/control paths are required. --original
verifies all arrays and non-runtime metrics against a same-source full run.
Algorithm changes instead use independent complete-workflow comparisons.

Bindings restore on exit. Resource interruption and numerical rejection are
distinct; neither authorizes an accepted rate or a relaxed check. No unit tests,
old source edits, state filtering or default solver adoption.
