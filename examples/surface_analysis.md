# Surface-layer accuracy analysis

Reads only complete configuration/source/hash-matched reference and dynamic
records. `include_extra=False` explicitly reads the original complete 20-case
protocol while a supplemental time-control job is still running; full assembly
requires the default, including that extra case. The comparison-pair list is
shared with the independently hashed integration cache. Comparisons retain
physical D and initialized state. The supplemental time-control configuration is
checked exactly through its separately hashed reproduction helper. Comparisons sample
fields at the same physical initial material coordinates. Gas velocity uses
initial volume weights (rho=1); fields use evolved test-volume weights. Changes
subtract each discretization's initial sample; total current and cumulative
local heat are also compared. Numerical field differences are not exact errors.

Coarse initial polygons lie within angular-refined polygons. The omitted finer
initial volume fraction is explicitly reported. Surface-vertex velocity/normal
velocity comparisons separately retain boundary errors that a bulk norm could
hide. Density/pressure/heat are the actual P0 cell fields, without smoothing.
Initial cell intersections form a common integration mesh, so discontinuous
cell fields and Q1-gradient traces do not cross integration regions. Positive
triangle rules integrate each overlap; the measured volume must match the whole
test polygon. Deformed volume weights include the actual Jacobian and r ratio.
Fourth/sixth-order quadrature checks assess remaining smooth integration error.
Nodal interpolation is never substituted for
exact flux freezing or a claim that a heated vacuum surface is resolved.
