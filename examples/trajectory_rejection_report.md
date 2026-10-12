# Refined trajectory rejection report

`run(localization, verification, destination)` renders the frozen density
localization and its complete byte-identical observation replay into PNG/PDF.
It requires the recorded primary17-step/step18 rejection and preserves every
input summary/raw hash with its own source archive. The separate expensive
failed-trajectory replay is explicitly not asserted complete by this report.

The four panels show all three prefixes' density error over the independent
linear signal, the final refined error partition across radial elements,
refined low/upper-quarter axial cell-mode fractions, and sampled internal
density jumps relative to the independent background. Quadrature20 is shown;
quadrature16/20 measurements remain in the underlying scientific records.

Fractional density jumps are multiplied by100 for the percent axis. Time is
displayed in microseconds. The zero initial interface jump is omitted from
the logarithmic jump panel only; all raw states and measurements are retained.
Nothing is fitted, smoothed or filtered. Initial error fractions have an
explicit precision/quadrature caveat, and spectra retain the aliasing caveat.

The report labels the old medium arithmetic difference and distinguishes
geometric localization from causal identification or a local conservation
budget. It is a research verification figure, not a formation or lab recipe.
Visual inspection of the generated PNG remains required after rendering.
