# Guard complete released-workspace refinement

Uses the unchanged guarded_tensor_refinement execution policy: F64/256 fluid,
selected field nr/nz, independent radial16,N128 reference, epsilon1e-4 and
original strict tolerance1e-10. Selects only the released-workspace full driver.
The control source archive explicitly includes this selection and its sources.

Own-worker peak RSS limit remains8GiB, free disk floor256MiB; native work that
does not yield to SIGINT exits only its own process after5s. No other process
is inspected or stopped. Original and replay use fresh directories; --original
requires all scientific arrays and non-runtime metrics to reproduce the supplied
same-variant full run. Different algorithms use the independent complete-workflow
comparison helper instead of claiming a same-source replay.

Bindings restore on exit. Resource and numerical failures remain distinct.
No equations, physical settings, bounds, budgets, source archives or user data
are modified. No unit tests are added or run.
