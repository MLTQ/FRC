# Refined amplitude and residual evidence report

`run(destination, audit_directory)` reads a terminal refined_control_audits
campaign and its unchanged physical campaign. It checks all complete replay
records, every recorded summary/raw/input hash and physical resource limits.
Each control is linked to its same-grid epsilon1e-4 baseline. All fifteen
scalar component norms, two vector norms, regional differences, independent
current/strong-curl comparisons and small-amplitude spectra stay in the JSON.
Twelve current error/reference squared integrals per accepted control must
agree exactly between the separate component and curl reconstruction audits.
Rejected tighter-bound protocols retain their failure message and bitwise
comparison of saved input array shapes, dtypes and bytes to the baseline,
without an error estimate.
This comparison establishes saved-input equality, not in-memory immutability
during a failed call.

The report preserves original signed residuals and instantaneous work closure
`abs(sum(P))/sum(abs(P))`, q16/q20 norm differences, discrete Parseval and curl
decomposition defects.
Region-wise quadrature differences stay separate: radius-band masks select
the original Gauss nodes and do not split cells at each band boundary.
Spectral fraction differences and current/strong-curl norm differences have
their own quadrature summaries. Spectral above-reference RMS is the component-relative
error times the square root of the fraction at cell modes abs(k)>64. Its small/
baseline ratio is descriptive; two amplitudes do not identify a noise source.
Cell modes fold aliases; inside-cell strong curl omits face distributions.

Standalone PNG/PDF plots each current component against amplitude, with hollow
strong-curl markers, and actual changes from the baseline. Tighter-bound
statuses appear separately, without plotting rejected rates. Inspect the PNG
after rendering. A fresh destination, immutable source archive, complete input
inventory and artifact hashes make the report reproducible and reviewable.

No plasma solve, modified state, interpolation of missing evidence, numerical
diffusion, filter, floor or relaxed tolerance is introduced. End-to-end here
means the complete signed-rate protocol, not time evolution or formation.
No unit tests. Prior source and calculation archives remain unchanged.
