# Species charge probe audit

Checks the frozen calculation closure and every referenced raw-state hash,
then replays orders 16 and 20 bitwise. Checks the independent continuum identity
between pressure work, entropy thermal work and enthalpy-gradient work. Reports
compatible-gradient annihilation and higher-quadrature agreement. Requires
potential-reconstruction errors to decrease on the measured refinement series.

The auditor records its own source hash alongside the already frozen probe
identity. This is a diagnostic audit, not a species integrator or full pressure
force-accuracy test. Fixed 1e-13 checks are scaled by the integrated absolute
thermal/pressure power, not by cancellation-prone nearly zero global work.
