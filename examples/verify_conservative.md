# Conservative-coupling report and plots

Builds a six-panel scientific figure and detailed Markdown report from the
recorded conservative measurement suite. Failed controls remain explicit and
are excluded from final-time convergence curves. Historical split shock curves
are labeled as previous measurements with their original fingerprint, not
presented as fresh identical-source runs. Whole-space response error is distinct
from exact fluid-cell conservation; source cases normalize it to actual work.
Limiter activity and physical versus numerical interface heat remain visible.
Fine-grid profile error being worse than the historical stabilized split result
is explicit, with both values and limiter activity taken from recorded data.
Boundary/time comparisons use separate readable percentage tables
for shock profiles and pulse observables.

`--reuse` accepts only source/environment-matched measurement cases. Rendering
source has a separate fingerprint, allowing report edits without repeating the
unchanged expensive solver measurements. Complete inputs and per-case progress
remain available if a longer verification run is interrupted.
