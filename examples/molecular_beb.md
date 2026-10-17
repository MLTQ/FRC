# Molecular total-ionization comparison

Run `python -m examples.molecular_beb FRESH_OUTPUT` from the repository root.
Load the exact attributed NIST inputs, calculate every independently tabulated
H2/N2/O2 cross section, and preserve the residuals without fitted corrections.
The reference prints 0.001 angstrom² increments, so its half-rounding interval
is explicit. Agreement verifies numerical reproduction of NIST's theory, not
an experimental validation or a complete discharge model.

For electron temperatures 1–100 eV, compare the full supported-energy Maxwellian
integral using Gauss–Legendre orders 8/16/32 and independent adaptive
Gauss–Kronrod integration. Each method resolves orbital thresholds and uses the
same source-defined cross section. The source formula is not independently
validated by comparing two quadratures; the external NIST tables serve that
purpose. The electron probability beyond 5000 eV remains in each integral record.

Outputs include immutable source/companion archive, exact input manifest,
tabular residuals and rate arrays, and scientific PNG/PDF. Terminal success
requires all published rounding intervals and a declared 1e-9 fine quadrature
comparison. Any failure is saved; parameters and tolerances are not tuned to
hide disagreement. This is a complete component calculation, with no unit suite
and no invented molecular reaction products or plasma trajectory.

The first `outputs/molecular_beb` run used the NIST introduction's rounded
explanatory constants. Five of 293 values exceeded the printed half-rounding
interval by up to 1.82e-5 angstrom²; that terminal failure is preserved. The
current source uses independently published NIST CODATA 2022 constants and
records them explicitly. The original tolerance and reference values are
unchanged; use a fresh output for the corrected calculation.
