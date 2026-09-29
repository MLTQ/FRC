# Moving warm species report

Runs the strict complete protocol/source/archive/history/physical-energy audit
before reporting. Produces separately labeled physical continuous and node/cell
error figures and species profiles against exact contact/simple-wave evolution.
Time-halving differences use the 2048-step same-grid trajectory, in Pa and
relative to uniform pe; ratios must be interpreted separately from spatial
continuum error. Half-CFL, reference-density-scale and radial-resolution
controls quantify their effects without changing base records.
The six stronger-gradient/common-time controls are audited as a separate full
collection with a checked driver hash. Both collections share the frozen
engine; report coverage is all 36 cases. Radial controls at a common time step
are separated from the base radial comparison, which also changes the CFL.
Convergence axes label only the five actual resolutions to avoid overlapping
automatic logarithmic tick labels in the exported figure.

Report, audit, physical-norm and independent magnetic-energy audit sources are
hashed separately from measurement sources, as are rendered figures. A completed
report verifies this bounded ideal gate only, not the full formation objective.
