# Configuration
`Material` specifies a fixed cylindrical conductor and reference plasma data.
Density and temperatures inform validity estimates; they are not evolved.
Ionization is metadata for warnings, not a kinetic source model. `Experiment`
validates the chamber, source, plasma and outer-domain relationships and supplies
a uniform integration time grid with actual dt <= requested dt. `output_stride`
controls field/topology snapshots; scalar histories retain every integration
step. A one-ULP tolerance prevents decimal ratios from adding an unintended step.
`persistence_threshold` defines the sweep's sampled persistent-flux class.
JSON round trips
explicitly reconstruct types; callbacks cannot be persisted. `default_coils`
is a small illustrative two-coil fixture, not a laboratory prescription.
