# Recorded conditional physical regime scan

Scans explicit electron density, temperature and physical length for assumed
Ar+, plus the original zero-current reference, proton comparison and partial
ionization warning case. No gas-pressure conversion, ionization equilibrium or
recommended laboratory setting is inferred. Gradient current B/(mu0 L) is an
explicit scale estimate, not computed formation current.

A separate prescribed vanishing-density profile reports local scale failure
as vacuum is approached at fixed temperatures. Its current is intentionally a
fixed scale to expose increasing required carrier drift; it is not a solved
edge or an input secretly inserted into moving dynamics. Exact zero-density
vacuum is outside the plasma diagnostic domain. Source hashes and all inputs
and warnings are saved in transport.json for reproducibility.
