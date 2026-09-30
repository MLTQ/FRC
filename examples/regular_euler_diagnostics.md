# Matched physical comparison and figure

Reads both frozen14-case mixed protocols. All state fields are independently
sampled on the same physical grid, with true species temperatures/pressures and
regular velocities. Successive mesh differences are convergence evidence,
not an exact unforced solution. Base time steps differ because the schemes have
different spectra; each has its own half/quarter-step controls.

For manufactured fields, the report separates exact nodal interpolation error
from evolution away from that interpolant. It retains fixed physical central
region errors. Between-point speed and grid/continuous energy discrepancy expose
the earlier representation problem instead of relying on conserved grid totals.
The six-panel Matplotlib figure uses actual saved data; no field is filtered or
renormalized for display. Source/audit hashes accompany the figure.
