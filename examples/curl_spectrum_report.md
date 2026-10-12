# Curl arithmetic and physical mode-spectrum figure

Loads terminal curl and cell-mode audits, checks spectrum array hashes, and
plots alternate arithmetic, complete radial-current error spectra, per-amplitude
component errors and exact-raw-curl comparisons. Repeated amplitudes from
different complete audits must have identical physical spectra and integrals.

Spectral energies for positive and negative cell modes are folded together
only for plotting. The underlying complete signed spectra remain unchanged.
Within/above-reference contributions are norms and combine in quadrature, not
by adding their heights. Input hashes and the figure/source hashes are recorded.
The source audit is still required to establish original physical fingerprints
and Parseval/quadrature/replay controls. Visually inspect the PNG after rendering.

The figure reports a filled-periodic verification fixture and unresolved
numerical accuracy, not a physical formation experiment or solver repair.
