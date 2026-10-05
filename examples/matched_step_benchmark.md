# Frozen nonlinear-step comparison

Runs explicitly requested complete steps from one signed independent continuum
mode, using the unchanged uniform-inertia plasma equations and either modal
Gauss solver. Freezes the entire imported source closure and companions, exact
initial/accepted states, actual stage changes, residual histories and runtime.
A rejected step records its type/message and last accepted physical state;
terminal never means all requested steps passed. Refuses output overwrite.

This is a performance/proposal control. It is not a physical convergence study
or an independently audited full mode trajectory. The default epsilon is 1e-4,
duration scale 100 microseconds and intervals 60, all explicit in the result.
No fixed-step retry, tolerance change or accepted-state repair is allowed.
