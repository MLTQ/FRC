# Reuse unchanged axial reconstruction inside streamed workflows

The first streamed interface audit calls magnetic_curl for each radial slice.
That frozen helper rebuilds the identical axial second-derivative basis on
every call. This wrapper constructs it once per complete interface audit and
reuses its array. Every radial basis and matrix expression is identical to the
original helper. The cache is reset between meshes and quadrature orders.

run dispatches the complete streamed_radial_rates workflow; all physical
solves, early raw records, volume/face diagnostics, norm reductions and strict
bounds remain. Source closure explicitly records this wrapper and both base
modules, including imports that use `from . import ...`. Temporary serial
bindings restore even on failure. Use separate processes for concurrent runs.

Compare full signed same-grid outputs against the frozen unbatched baseline
and the first completed streamed baseline. Exact states/rates/face arrays and
solver diagnostics must agree bitwise. Diagnostic differences remain explicit;
the cached and uncached same-batch audits should also match every fingerprint
and metric. Verify new workflows through full fresh replays. No unit tests.

This is diagnostic reuse, not a new solver, preconditioner, spatial method,
physical approximation, boundary change or numerical diffusion.
