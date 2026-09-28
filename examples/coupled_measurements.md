# Fixed-chamber verification measurements

Collects force, passive advection, Gaussian diffusion, fixed-grid temporal order,
infinite-space vacuum comparison, and prescribed two-coil pulse data. Every
dynamic row records full mesh, gas parameters, bias, coils and serialized waveform,
plus energy/current/divergence diagnostics. Initial gas is rho=1 kg/m³, p=1 Pa;
only the advection fixture has initial vz=0.5 m/s. Mathematical initial fields
and comparison windows are specified by the benchmark module and report.

`pulse_measurement` samples 21 uniformly spaced times through the complete 0.1 s
pulse and optionally stores final fields on their distinct fluid/magnetic grids.
Snapshot time boundaries are included in the accepted step schedule. The pulse
is examined under separate space, outer-domain and CFL refinements; no topology
or formation criterion is inferred from these integration checks.

`all_measurements` runs the bounded study, printing actual completed cases so a
long run's progress is visible. Failures propagate; partial trajectories are not
silently reported as completed comparisons. Presentation and provenance live in
`verify_coupled`, keeping measurement code independent of plot style.
