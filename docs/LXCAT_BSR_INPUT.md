# Private LXCat BSR argon input

On 2026-10-16 the user explicitly authorized “Accept for local research.”
The LXCat terms were accepted through the official browser workflow, and the
BSR electron–argon archive was downloaded successfully. This supersedes older
“permission pending” notes; it does not authorize redistribution or commercial use.

Private storage, outside the repository:
`/Users/max/Library/Application Support/FRC/personal-atomic-inputs/lxcat-bsr-argon-2026-10-16/`

The original archive and unchanged text/XML exports are preserved alongside
`provenance.json` and a screenshot of the accepted download workflow. The text
contains 33 process blocks: one elastic momentum-transfer, 31 excitation and one
ionization. The header identifies BSR-500. The ground-neutral ionization channel is now an explicit optional input to
the homogeneous atomic source; the completed comparison below measures its
effect. Excitation state mapping and suitability for other reactions still
require review. Do not infer excited-parent or final-state branching data
from the number of excitation channels. The moving formation solver and visualization remain unchanged. The user
authorized this bounded continuation after the earlier pause.

## Attribution and restrictions

Reference: **BSR database, www.lxcat.net, retrieved on October 16, 2026.**
The provider header additionally requests:

- O. Zatsarinny (2006), “BSR: B-spline atomic R-matrix codes,” Computer Physics
  Communications 174, 273.
- O. Zatsarinny and K. Bartschat (2013), “The B-spline R-matrix method for atomic
  processes: application to atomic structure, electron collisions and
  photoionization,” Journal of Physics B 46, 112001.

Preserve applicable process-specific references when using individual channels.
[LXCat terms](https://us.lxcat.net/instructions/how_reference.php) require database,
site, retrieval date and original-source attribution, and prohibit third-party
distribution and commercial use. Raw tables and reconstructed provider rate
arrays must not be bundled in Git, source archives or shared artifacts.
[Database](https://www.lxcat.net/BSR),
[official download workflow](https://us.lxcat.net/data/download_set.php).

## Exact input identities

| Private file | Bytes | SHA-256 |
|---|---:|---|
| download.zip | 68070 | `0d7058707eeb0ac25ee7d537d6e119bb8a42c9d586a689cb2dc6233e5b0c6123` |
| cross_sections.txt | 183265 | `ed851c2f5f9292ff46ff890069a28b3f4c4edba1068ae130c0a3dd6080758f44` |
| cross_sections.xml | 174279 | `c7ecf9762ccedac10113cc77fddbbf7209e8d4ad2c97cfa4e38241a3b6b4ea3b` |

On resumption, inspect the actual target/product identities and cross-section
energy domains; map only defensible channels into the existing atomic source,
retain explicit data gaps, and compare complete trajectories with independent
baselines. Track this through frc-transport.5.3 and frc-transport.5.4. Do not
start new simulation work solely because this access permission was granted.

## Ground-neutral ionization comparison completed (2026-10-16)

The user asked what LXCat would do and authorized further work. New generic
`lxcat_cross_sections` integrates the explicitly supplied cross sections for
zero-drift Maxwellian electrons. `cross_section_ionization` replaces only
Ar(ground) ionization in the existing `AtomicNetwork`; it does not double count
the old fit, alter physical atomic energies or assume inverse recombination.
All other atomic channels stay unchanged. No raw tables/rate arrays are bundled.

Run73644 EXIT0, all12 full3.5us trajectories in101.749s, source
`1fe8ad2f3d218e6f32ddc32ebe03ff96aa850b9bf85ad2ac1f0397b39e5d087c`.
Outputs: `outputs/bsr_atomic_source`, all accepted states, source archive,
summary and inspected scientific PNG/PDF. BSR runs use350/700/1400steps;
each has a matched old-rate control at1400steps. These are homogeneous source
histories with assumed seeds/heat, not physical gas preparation or formation.

At10eV, the BSR ground-ionization coefficient is0.5930473 times the former
Voronov fit. The effects at the finest2.5ns step are:

| Initial condition | Final electrons/nucleus: old → BSR | Final Te: old → BSR [eV] | Peak Te: old → BSR [eV] |
|---|---:|---:|---:|
| Initially singly ionized | 1.04927679 → 1.04927594 | 2.31777975 → 2.31778354 | 10 → 10 |
| 10% seed, no heating | 0.13564697 → 0.12490665 | 2.38680445 → 2.62839189 | 10 → 10 |
| 10% seed,20eV/nucleus heat | 0.68390234 → 0.62667761 | 2.53628064 → 2.71224324 | 17.1787421 → 21.1041846 |

The heated case's maximum between-data-source temperature/charge differences
are3.99324eV /0.0572247 electrons per nucleus. BSR5→2.5ns differences are only
0.00112014eV /8.96831e-6;10→5ns errors are aboutfour times larger. Across all
cases, number/charge discrepancies are<=2.72e-13 /1.61e-13 and energy balance
error<=3.78e-12eV/nucleus. The independent prescribed-heat quadrature discrepancy
is separately<=0.00052366eV/nucleus; it is not the same conservation diagnostic.

The independent finite-domain integral for sigma(E)=a+bE has maximum relative
quadrature errors7.23e-5,1.79e-8,2.22e-16 at4/8/16points per interval. Actual
BSR8/16point differences are<=1.06e-13 over the sampled1–40eV temperatures.
This verifies integration of a linearly interpolated table, not the accuracy
of the underlying table or the Maxwellian/zero-drift assumption.

The ionization table spans15.76–300eV. No high-energy extrapolation is invented.
Its threshold differs from the independent physical15.7596119eV energy by
0.0003881eV, retained explicitly. The heated fine history's maximum Maxwellian
speed-weighted fraction above300eV is1.02023e-5. That is NOT a bound on missing
reaction rate. No validated Maxwellian T interval is supplied: the replacement
is marked unspecified/outside-fit-range in existing diagnostics. Ground-final
branching, excited-parent omissions and produced-photon/escape limitations remain.

## Why the excitation table was not silently merged

All31 excitation records have ground-neutral targets; they do not fill missing
excited-parent ionization. For the named four4s states, BSR/NIST energy differences
are<=0.000646eV. For the ten named4p states, differences are0.1121–0.1670eV even
though statistical weights match. This is an energy-convention/atomic-data issue
requiring explicit treatment; neither the table axes nor physical NIST levels
were shifted to conceal it. Higher-level and lumped `nl>5s` records also require
consistent level/radiative modeling before integration.

Next: compare the competing ground-ionization sources with independent physical
reference data and resolve excitation identities/energy conventions before
choosing a default. Keep the existing moving-boundary taskfrc-transport.5.5.1
and full formation/transport gates open. This completed taskfrc-transport.5.6
is an opt-in source comparison; it does not resolve the free-boundary failures.
