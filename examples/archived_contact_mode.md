# Explicit independent contact reference selection

Selects a hash-checked entry from a complete independent continuum spectrum
archive. Retains its exact radial wavenumber, axial coefficients, growth,
physical scale and original real-Nyquist convention. The inherited analytic
J0/J1 reconstruction, original full electric field and current are unchanged.
Canonical real sign is the same maximum-component convention as the old first
radial reference. No production eigenvector enters this reference.

`ArchivedContactMode(directory,index)` defaults to refined radial16,N128 entry1.
The original radial16,N64 reference retains a measurable real-Nyquist curl
defect; it is available only by explicit archive selection for comparison.
The source summary must have a `continuum` list with the
same explicit raw-array schema as stationary_contact_audit. Records both raw
and summary hashes. This opt-in class leaves all old reference sources frozen.
