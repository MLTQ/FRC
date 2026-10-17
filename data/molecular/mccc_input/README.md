# Private MCCC molecular inputs

On2026-10-17 the H2 and D2 ground-electronic-state collections were downloaded
from the [authors' MCCC database](https://mccc-db.org/) for this local research
project. The authors prohibit third-party redistribution and commercial use.
Raw ZIPs remain outside Git in `/private/tmp/frc-mccc-personal/`; preserve them.
Only provenance, coverage and interpretation notes are recorded here. Cite the
dataset publication and database/retrieval date when using these data.

`manifest.json` records exact archive hashes, sizes, request fields, local paths,
channel-family counts, initial-level coverage and numerical-file observations.
The stable download operation is POST to `https://mccc-db.org/cgi-bin/archive.php`
with its recorded `fileName` and `zipName`; the response supplies a temporary
ZIP path relative to the site. No account, agreement-acceptance control or
outgoing private research data was involved. The ZIPs were read directly and
not extracted or modified.

| Ground-state archive | Raw cross-section members | Initial vibrational levels | Ionization coverage |
| --- | ---: | --- | --- |
| H2 |6316|v=0..14; rotational collection only vi=0|v=0..14|
| D2 |12385|v=0..21|v=0..20; v=21 missing|

Counts include two empty H2 rotational files; they are not usable tables.
Electronic-excitation collections contain both individual final levels and
aggregate bound/continuum/total curves. These counts must not be interpreted
as disjoint reactions. The D2 archive has no direct X-to-X excitation or
rovibrational collection; its ERD curves are a different physical process.

The H2 Ni=30→32 and31→33 files are empty. D2 b3Su dissociative excitation
from vi=15 has commas after the first four incident-energy values, followed
by whitespace-only rows. The original is retained; a future reader must handle
that documented delimiter explicitly without changing numerical values. Other
nonempty tables parse as two finite columns with increasing incident energy and
nonnegative cross section. Integrity checks do not establish physical accuracy.

## Physical interpretation before coupling

[Electronic excitation](https://doi.org/10.1016/j.adt.2020.101403) supplies
resolved bound final states and separate dissociative contributions, with
isotope-specific vibrational energies. Use those identities and energies; do
not make D2 a mass-scaled H2 alias. Final table energies differ by channel:
the D2 collection ends at150–500eV. No extension outside a table is implied.
H2 also includes low-energy X-to-X data identified by the current website as
submitted to Physical Review A(2026); its publication remains to be established.

[Total ionization](https://doi.org/10.1016/j.adt.2023.101573) combines
dissociative and nondissociative outcomes. It cannot all create an intact
molecular ion or specify electron/fragment energy sharing. The D2 vi=21 omission
agrees with the paper's21/22 ground-level coverage. Tables end at1000eV.
Fitted thresholds can denote substantial onset rather than the physical lowest
product threshold; do not use fit parameters as energy-level definitions.

[ERD](https://doi.org/10.1016/j.adt.2022.101534) combines electronic excitation
and subsequent radiative cascades. It is not direct vibrational excitation and
must not be added to an explicit excitation/cascade network as an independent
channel. Ground-state energy difference alone omits photon/electron energy
partition. D2 radiative dissociation is not supplied by these bound-final ERD
tables. Their common upper incident energy is500eV.

Next construct an attributed state/energy network and select disjoint processes,
with explicit photon and fragment reservoirs. Resolve direct D2 vibrational
collisions, missing ionization levels, product branching and angular/energy
moments before moving coupling. No H2/D2 discharge model is enabled by this
acquisition. Independent physical comparisons and complete coupled runs remain.
