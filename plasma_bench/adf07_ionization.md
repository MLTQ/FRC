# Explicit local state-resolved ionization input

Reads ADF07 reaction blocks from an explicit local file, retaining source
initial/final metastable identifiers and effective ionization potentials.
There is no downloader or embedded provider data. Input path/hash are recorded;
the raw personal OPEN-ADAS file remains outside the repository and archives.

The [official format](https://open.adas.ac.uk/man/appxa-07.pdf) specifies eV
temperatures, cm^-1 effective potentials and cm^3/s rates. The reader converts
only rates to SI and accepts ordinary E or Fortran D exponents. Positive rates,
increasing temperatures, unique reactions and adjacent charges are required.
Unsupported records reject. `coefficient` interpolates log(rate) in log(T)
within the supplied interval; it does not extrapolate or invent zero rates.

Source metastable numbers are not level-array indices. Mapping requires an
explicit physical state identification by the caller. In particular final
metastable code 0 is unresolved; it is not a supplied ground-state branching
fraction. Any model collapsing it to ground must disclose that assumption and
the unavailable final excitation energy. Source effective potential is retained
separately from the independent energies in the model's energy inventory.

`atomic_network` accepts explicit (charge, source-level, data-row) mappings;
the argon example identifies the three initial LS terms from source comments
and the excitation input's configurations, term multiplicities and weights.
