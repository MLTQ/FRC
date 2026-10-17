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
HTML/ASCII cross-section download links returned HTTP500 during this review;
no independently tabulated cross sections were imported or compared.

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
before enabling a gas as a lab-facing prediction. None of these new molecular
gases is currently enabled by this source review.
