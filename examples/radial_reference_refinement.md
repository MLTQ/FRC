# Independent axial reference refinement at fixed radial mode

Builds the original continuum Fourier/Bessel matrix at radial16 and axial
64/128/256 points. Records the full operator, spectrum, real leading vector,
physical scales, divergence and eigensolve backward error. No production fluid,
field, current, projection or time-stepping operator is imported.

The radial wavenumber stays fixed. Axial eigenvalue and physical eigenfunction
convergence must be compared separately; a leading eigenvalue alone cannot
validate a different numerical mode. Output retains the explicit continuum
archive schema for ArchivedContactMode, source/runtime and raw hashes. Fresh
complete replay is required. Existing independent spectra remain unchanged.
