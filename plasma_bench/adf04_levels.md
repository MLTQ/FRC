# Explicit personal-input excitation reader

Reads the untagged type-3 Maxwellian excitation subset of ADF04 from an explicit
local path. It performs no download and embeds no provider data in the code.
The raw input hash/path are retained; callers must respect the input's own terms.
The inspected OPEN-ADAS inputs are kept outside the repository for personal local
research, and are not copied into source or experiment archives.

The official [ADF04 specification](https://open.adas.ac.uk/man/appxa-04.pdf)
defines shortened Fortran exponents, kelvin temperatures, effective collision
strengths, A coefficients and optional Born-limit columns. Unsupported process
format types reject. Additional tagged processes are returned verbatim with
source line numbers; the population component must explicitly declare any
excluded process, rather than silently dropping physics. Literal
zero/tiny source values are retained; an absent row is not a demonstrated zero.
No high-temperature limit or extrapolation is manufactured from the optional
Born column.

Source xj is (statistical weight-1)/2 in LS coupling and J in intermediate
coupling. Thus g=2*xj+1 in either case, but the reader does not label every xj a
physical angular momentum. Configurations, multiplicities, orbital L and source
energies remain explicit. IDs must be contiguous and transition pairs unique.
`read_levels` returns metadata and arrays for `collisional_levels`; it does not
infer a particular coupling scheme, add unlisted transitions or fit atomic rates.
