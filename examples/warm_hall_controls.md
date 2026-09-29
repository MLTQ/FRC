# Follow-up warm accuracy controls

Six explicit supplemental configurations retain the main 58-case protocol.
Pressure-only nr64 uses 160/320/640 steps and nr128 uses 320 steps to separate
the nonmonotonic maximum local ion-pressure error from time and mesh error.
The slow branch at nr16 repeats with half fluid CFL and one tenth wave amplitude
to distinguish time splitting and nonlinear effects from spatial errors.

Copies of complete base configurations change only named settings. The control
driver has its own frozen hash, checked after each evolution and at assembly,
in addition to the unchanged transitive engine identity. Independent pressure/
slow groups write separate files and require all six records for assembly.
No case or failure is silently omitted. Higher resolution does not itself
establish acceptable accuracy; all original coarse data remain in the report.

Observed result: the nr128 pressure case raises before archival; its group has
three saved successful controls. The slow group saves two. Dedicated
`warm_pressure_failure` captures the original failed case and two finer-time
rejections, and `verify_warm_hall_controls` checks complete success/failure
coverage. `--assemble` consequently rejects this incomplete-success collection;
do not create a misleading all-success `measurements.json` for these controls.
