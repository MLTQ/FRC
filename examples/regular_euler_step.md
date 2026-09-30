# Independent forcing/work for active regular-coordinate DOFs

The analytic source is the original regular Euler residual. Inactive axis u/
Omega slots receive no forcing. Original physical energy is differentiated
independently by complex step at the **actual** RK stage state to measure source
work. Integrals use actual angular momentum2x*l, not the unweighted l sum.
Diagnostic work/source ledgers never alter an evolved state. Unforced steps
call the real regular driver unchanged; every final forced state is validated.
