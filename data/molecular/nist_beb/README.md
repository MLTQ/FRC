# NIST molecular total-ionization inputs

Original public NIST ASCII tables for H2, N2 and O2 are preserved byte for byte.
`manifest.json` records URLs, retrieval time, units and SHA-256 hashes. The
orbital parameters are inputs; the separately generated cross-section tables
are an external numerical reference for the BEB implementation. Cross sections
are in square angstroms (1 angstrom² = 1e-20 m²), energies in eV.

The legacy CGI links use a bare query, such as `bebcsdwnload_ascii?H2`. Appending
an equals sign changes the request and caused the earlier reader's HTTP 500
errors. The actual site links were verified in the browser before acquisition.

Sources: [NIST explanation](https://physics.nist.gov/PhysRefData/Ionization/intro.html),
[equations](https://physics.nist.gov/PhysRefData/Ionization/Eqs/latex.html),
[H2](https://physics.nist.gov/cgi-bin/Ionization/table.pl?ionization=H2),
[N2](https://physics.nist.gov/cgi-bin/Ionization/table.pl?ionization=N2),
[O2](https://physics.nist.gov/cgi-bin/Ionization/table.pl?ionization=O2).

These are theoretical total-ionization data, not a complete molecular collision
set or a new experimental calibration. The supplied orbitals all have Q=1,
no double-ionization multiplier and no special denominator. Only those declared
channels are implemented. No D2 data or isotope substitution is supplied.
Product identities/energies, dissociation, attachment, internal states and
transport must be sourced before coupling these totals to a molecular plasma.
