# Available gases and required physical models

User scope,2026-10-16: H2,D2,N2,O2 and common available gases; retain argon as
the current numerical baseline. Provisional chamber/fill inputs are roughly 1 m
length, 15 cm diameter and 0.1–100 mTorr. No operating recommendation follows from
the present planar argon component calculations. Bead `frc-transport.6` stays open.

The older `material_experiment.GASES` table already accepts Ar, He, Ne and H2
for an assumed fully singly ionized inventory. Its H2 entry assigns two protons
per original molecule, assumes complete dissociation, and evolves neither
molecules nor preparation energy. This is an initialization assumption, not an
implemented hydrogen discharge model. D2, N2 and O2 are not aliases in that table.
The current planar reacting calculation uses the argon ground-state ladder;
it is a different component from that older prescribed-coil experiment.

## Primary data sources reviewed

[NIST's electron-impact database](https://pml.nist.gov/PhysRefData/Ionization/intro.html)
provides the Binary-Encounter-Bethe total-ionization model and molecular orbital
inputs. Its pages for [H2](https://physics.nist.gov/cgi-bin/Ionization/table.pl?ionization=H2),
[N2](https://physics.nist.gov/cgi-bin/Ionization/table.pl?ionization=N2) and
[O2](https://physics.nist.gov/cgi-bin/Ionization/table.pl?ionization=O2) are accessible.
Total ionization alone does not specify product branching, fragment energy or
the other collisions needed for a discharge. The database distinguishes total
and differential information and declares the model nonrelativistic. Its
HTML/ASCII links initially failed because the reader changed their bare query
strings. Browser verification and direct downloads resolved this: the original
H2/N2/O2 orbital and cross-section tables are now preserved with URLs and hashes
in `data/molecular/nist_beb`. They provide an independent numerical reference
for the BEB implementation, rather than an experimental calibration.

The [MCCC database](https://mccc-db.org/) exposes isotope and vibronic-state
selection. The authors' [IAEA contribution](https://conferences.iaea.org/event/359/contributions/29407/)
describes vibrationally resolved hydrogen/isotopologue calculations. This is a
candidate source for distinct H2/D2 channels. MCCC restricts redistribution and
commercial use and requests dataset-publication plus database citations. No
MCCC data were downloaded or bundled here; review each selected dataset's
coverage and applicable terms before local use.

The [LXCat authors' N2/O2 comparison](https://us.lxcat.net/notes/index.php?download=lxcat_team10)
distinguishes swarm-consistent compilations from partial recommended cross
sections. It includes rotational/vibrational/electronic excitation, dissociation,
ionization and oxygen attachment, with dataset-specific restrictions. Thus a
named database alone is insufficient provenance: record version, process
meaning, threshold, units, electron distribution, validity and primary source.

## Implementation requirements, not yet supported features

| Gas | Work required beyond the current argon source |
|---|---|
| Ar | Finish local momentum/thermal transfer and physical transport/formation convergence; excited-state channels remain incomplete. |
| H2 | Molecules, atoms, molecular/atomic ions and applicable negative carriers; rotational/vibrational energy, dissociation, recombination and product branching. |
| D2 | Isotope-specific masses, level energies and reaction data; do not relabel H2 or scale only the mass. |
| N2 | Molecular internal-state energy, dissociation products, applicable ion chemistry and attributed transport data. |
| O2 | Molecular/atomic species plus applicable negative ions, attachment/detachment and dissociation/recombination energy. |

These are project model-design requirements, not a claim that the listed sources
supply a complete reaction set. Channel selection must follow the simulated
distribution and regime. A total ionization curve must not silently be assigned
to one molecular-ion product. Chemistry must conserve nuclei, charge, momentum
and declared energy reservoirs with physical masses and explicit radiation.
Use independent experimental/swarm comparisons and complete coupled refinement
before enabling a gas as a lab-facing prediction. None of these molecular gases
is enabled as a coupled discharge model by the component work below.

## Implemented total-ionization component

`plasma_bench/molecular_beb.py` loads the attributed H2/N2/O2 orbital constants
and implements NIST's declared ground-target BEB formula. It accepts only the
supplied unmodified Q=1 orbitals, with no double-ionization multiplier or
invented product branching. D2 is explicitly unsupported by this dataset.
Energy inputs are limited to the 0–5000 eV comparison range. The optional
Maxwellian integral assumes stationary targets and zero electron drift, returns
only the supported-energy contribution, and reports the omitted electron
probability without renormalizing the distribution.

The complete `examples.molecular_beb` calculation in
`outputs/molecular_beb_codata` reproduces all 293 independent NIST tabulated
values within their printed 0.001 angstrom² rounding. Maximum absolute
differences are 0.0004982 / 0.0004983 / 0.0004889 angstrom² for H2/N2/O2.
Across electron temperatures 1–100 eV, maximum relative differences from
independent adaptive quadrature fall from 1.20e-6 at order 8 to 1.14e-11 at
order 16 and 5.58e-16 at order 32. The largest electron probability above
5000 eV is 1.55e-21; it is not a bound on the unknown rate outside that range.
Raw numerical comparisons and the inspected scientific PNG/PDF are retained.

The first calculation using the introduction's rounded explanatory constants
is preserved as a terminal failure in `outputs/molecular_beb`. The corrected
version uses independently published NIST CODATA 2022 Bohr radius and Rydberg
energy, keeping the original data and tolerance. Source identity:
`ef2c7cb935d2cd0b39f52b902ab554dbd0beffee77baacbc5cb6984c29e9b615`.
Bead `frc-transport.6.1` covers this component. Molecular products, energy
partition, excitation/dissociation/attachment, D2 data and moving coupling stay
open under `frc-transport.6`.
