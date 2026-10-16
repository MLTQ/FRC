# Personal LXCat cross sections and Maxwellian integrals

Reads explicit user-supplied LXCat text files with eV/m² two-column tables.
No downloader, provider tables or reconstructed rate arrays are included. Keeps
literal process names, parameters, comments, energy grids and input SHA-256.
Supported process records are elastic, effective, excitation, ionization and
attachment. Rejects malformed units, nonfinite/negative values and unordered
grids. Reading an excitation record does not map it to a modeled atomic level.

For isotropic zero-drift Maxwellian electrons, T and E in eV, the finite-domain
rate in m³/s is

    k(T) = sqrt(8 e/(pi m_e)) T^(-3/2)
           integral[E_min,E_max] sigma(E) E exp(-E/T) dE.

`MaxwellianIntegral` linearly interpolates each source interval, then uses
positive Gauss–Legendre quadrature within that interval (default 16 points).
It neither shifts thresholds nor extrapolates below/above supplied energies.
The result is explicitly a truncated rate, not a claim that unsupplied cross
sections are zero. A coarse source grid can retain physical interpolation error
even when integration-order differences are small.

`electron_flux_fraction_above_table` returns `(1+x)*exp(-x)`, x=E_max/T.
This is the Maxwellian speed-weighted electron fraction above the table,
NOT a bound on missing reaction rate: the unsupplied cross section is unknown.
No fit validity interval is inferred from this number. Finite relative electron/
neutral drift or non-Maxwellian distributions require a different integral.

The BSR source is used only by explicit opt-in callers. Numerical diffusion in
this component is absent; its approximations are table interpolation, finite
energy support, Maxwellian averaging and quadrature. The companion example
compares quadrature to an independent closed-form integral and runs the existing
coupled atomic evolution. No unit suite or new reporting infrastructure.
