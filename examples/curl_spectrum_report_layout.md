# Clear legends in the archived scientific report

Reuses the unchanged full data/hash/plot workflow. Immediately before closing
the figure, extends panelC/D vertical limits and places the panelD legend above
its bars. Values, spectra, labels and physical conclusions remain unchanged.
The original figure/source archive is preserved in a separate directory.

Temporary serial close/provenance bindings restore on success or failure. The
new immutable archive includes both plotting sources, and the original report
hashes the final saved PNG. This module is for serial report generation only,
not concurrent plotting. Inspect the final PNG; no physics rerun is required.
