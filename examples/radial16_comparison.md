# Frozen radial16 controls

Requires the already passing old/new baseline to match every saved array
bitwise, including full RHS and physical-face integrals. For same fine grid at
amplitudes1e-4/1e-5, saves centered full-RHS coefficient differences and reports
Euclidean relative coefficient norms, explicitly distinct from physical norms.
Physical quadrature16/20 convergence compares each relative volume, axis and
wall norm in all complete cases; max errors sampled at different points are
not treated as quadrature integrals. Failed cases remain excluded and failed.

Checks every input raw hash and terminal batch before freezing this comparison.
A full fresh rerun through wall_diagnostic_replay is required. These controls do
not close the physical accuracy or full evolving-response gate.
