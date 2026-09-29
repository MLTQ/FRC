# Scientific resistive figures

`error_figure` distinguishes continuum field/current/heat error, surface/axis
error, independent stiff-mode time control, actual resistive relaxation, changed
physical diffusivity and moving source-work time refinement. Explicit resolution
and time ticks avoid logarithmic label collisions. Different D changes physics.

`field_figure` reuses the unchanged material-geometry/flux display helper for
actual saved pressure and volume-averaged toroidal current. It colors no exterior
fluid and labels the source fixture as initialized, not formed. Heat histories
and velocity projected on the initial surface radius complement field plots.
The current is a quadrature average, not a recovered pointwise peak. Both PNGs
require visual inspection, and all reused plot/geometry helper hashes belong
in the report provenance without changing historical study fingerprints.

Both completed-study PNGs were visually inspected at the resistive checkpoint:
labels, units, legends and ticks are readable; actual current/pressure and the
large surface-flow mesh differences remain visible. No cropped panels or
overlapping tick labels were found.
