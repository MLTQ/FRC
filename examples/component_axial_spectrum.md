# Physical component error by axial cell mode

Reconstructs both complete signed physical responses and the independent
continuum reference at every original quadrature16/20 point. All original
field/reference/weight fingerprints must match. Each scalar component error,
amplitude difference and reference is reshaped into (radial, axial cell,
intra-cell Gauss node) and Fourier transformed only along the cell index.

For n axial cells, fft/n is used and spectral squared magnitudes are multiplied
by n and the original radial/intra-cell weights. Parseval equates their sum to
the full original physical squared integral; defects are explicitly measured.
No mode is discarded or changed in a state or field. Signed mode indices,
complete spectra and grouped fractions are retained for every component.

These are translation-by-one-cell modes, not unaliased continuous Fourier
harmonics: harmonics separated by n fold into the same cell index while the
intra-cell coordinate is retained. The above-reference group has abs(k) >
Nreference/2. Groups overlap (e.g. upper-quarter lies within above-reference)
and must not be added together. High-mode error is diagnostic evidence, not
proof of a cause or authorization for filtering. Quadrature convergence and
complete fresh replay remain necessary. No unit tests are used.

The same frozen source, accepted-input hashes and actual positive/negative
responses are required. Pre-existing output directories reject; all arrays
and non-runtime metrics must independently replay before adoption of results.
