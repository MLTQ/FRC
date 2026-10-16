# Argon direct-ionization inputs

`cfit.dat` and its field description `cfit.txt` were downloaded unchanged from
[D. Verner's atomic-data server](https://www.pa.uky.edu/~verner/col.html),
which distributes the Voronov (1997), ADNDT65,1 rate fits. The complete source
table is retained; `argon.json` selects its18 argon transitions. The rate is
A*(1+P*sqrt(U))/(X+U)*U^K*exp(-U), U=threshold_fit_ev/Te_ev. A is tabulated in
cm³/s and must be multiplied by1e-6 for SI. Tmin is eV, Tmax is keV in the
original table. This implementation retains the analytic exponential tail,
not the optional U>80 early return in the historical Fortran convenience code.

`nist_argon.html` is the unchanged NIST ASD ionization-energy query response
recorded in `argon.json`. Its18 successive ground-state energies, displayed
uncertainties, bracket/parenthesis annotations and references are retained.
Number formatting spaces are removed only when converting the displayed
energy to a number (including the final4120.66559 and4426.22407eV entries).
The model's atomic energy uses these physical energies, separately from the
rate-fit threshold parameters. Not all NIST entries are direct measurements;
the original annotations and citations distinguish evaluated/estimated values.

In particular, the archival Ar3+ fit uses52.3eV while NIST lists59.58eV.
The fit is not silently retuned; this data discrepancy must remain visible in
regime diagnostics before physical adoption. The public table is an initial
direct-ionization source, not a modern complete argon discharge data set.
Excited states, metastables, excitation-autoionization, recombination,
radiation, opacity and non-Maxwellian electron distributions require further
physics/data. No OPEN-ADAS restricted data is incorporated here.

Each raw source hash and retrieval date is in `argon.json`. These local
research inputs are archived for reproducibility; no publication is requested.
