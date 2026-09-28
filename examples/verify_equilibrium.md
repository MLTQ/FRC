# Closed-flux equilibrium experiment report

`run_case` measures the independent continuum equilibrium's initial force error,
finite-time velocity/flux drift, mass/energy and divergence. It catches only
explicit nonphysical-state errors and records the rejection, message and actual
time reached. Unexpected software errors still fail the run. Failed partial
trajectories are never plotted or compared as completed equilibrium evolutions.

`verify` runs two declared background pressures, multiple spatial resolutions
and time refinements of a resolved 48-cell equilibrium. It preserves physical and
numerical parameters, source fingerprint, versions and field arrays. `plot`
shows initial topology and artificial flow; `report` records limitations without
claiming formation, stability, confinement or arbitrary equilibrium preservation.
The cited spectral-element reference motivates the equilibrium class; it is
not claimed as the source of this independently derived periodic series.
`run_case(positivity=False)` retains the formerly failing reconstruction as an
explicit reference. The report compares it with enabled limiting, records
accepted-stage activation fractions and never treats a completed limited run
as proof of arbitrary full-update positivity or adequate spatial resolution.
