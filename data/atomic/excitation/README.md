# Low-lying neutral argon radiation data

Retrieved 2026-10-16 from the public NIST Atomic Spectra Database (ASD 5.12,
2024), Kramida, Ralchenko, Reader and NIST ASD Team,
https://doi.org/10.18434/T4W30F. Query at
https://physics.nist.gov/PhysRefData/ASD/lines_form.html:
Ar I, upper energy <=109000 cm^-1, only lines with transition probabilities,
all-vacuum wavelengths, cm^-1 energies, configurations/terms/J/g/references,
CSV output. The original unmodified response is `nist_ar1_low_lines.csv`.

`argon_low_radiative.json` is a lossless numerical transcription of the selected
level energies, degeneracies and 31 transition probabilities, with original
accuracy codes and reference keys. The 15 unique energies are sorted; Paschen
labels are ground, 1s5/1s4/1s3/1s2, 2p10 through 2p1. Configuration, term and J
are retained so names alone cannot silently determine level identity. CSV
spreadsheet quoting is removed as text, never executed. Source SHA is recorded.

NIST identifies primary transition compilations Fuhr & Wiese (1996),
Wiese et al. (1989), and Morton (2003); rows retain T5223, T1242n, T7574.
The two resonance rates are 5.32e8 and 1.320e8 s^-1 (104.822 and106.666 nm).
Do not replace current ASD values with older Handbook values: some visible-line
probabilities differ. The model computes vacuum frequency from level energies
for consistent energy bookkeeping, retaining tabulated wavelengths separately.

This is a bounded set, not a complete argon atomic model. This directory contains no collision
cross sections, higher-state cascades, Ar II excitation, metastable destruction,
excited-state ionization, pressure broadening or line-transfer data here.
No missing transition is certified zero. Natural plus thermal widths alone
cannot establish an experiment's total line shape.

BSR electron-impact cross sections were acquired privately from LXCat after
the user authorized local research on 2026-10-16. They remain outside Git; see
`docs/LXCAT_BSR_INPUT.md` for hashes, attribution, restrictions and next-use limits. Never treat other datasets' non-Maxwellian E/N-dependent rate tables
as Maxwellian temperature rates. CHIANTI's public Ar list starts at Ar III and
does not supply the missing neutral/Ar II rates.

## Personal local excitation inputs (2026-10-16)

Public OPEN-ADAS inputs `arlike_mom08#ar0.dat` (15 Ar I levels) and
`cllike_dcg07#ar1.dat` (40 Ar II LS terms) were used for local personal research.
They remain outside the repository in `/private/tmp/frc-openadas-personal/`;
raw files and reconstructed rate arrays are not bundled with code or archives.
Exact links, identities, assumptions and completed results are documented in
`docs/COLLISIONAL_LEVELS.md`. The generic reader requires explicit local paths.
These are supplied Maxwellian effective collision strengths, separate from the
subsequently authorized private LXCat BSR cross-section download.

The Ar I excitation study explicitly uses this directory's current NIST
radiative data with matched ADAS collision levels; it excludes the source's S
ionization row only for the fixed-charge component. Ar II retains supplied
term-averaged energies/A values and is not fine-structure line transfer.
This progresses excitation physics but is not a complete cooling/charge model.

The follow-on joint source also accepts personal ADF07 `szd97#ar_ar1.dat` for
Ar II 2P/4D/4F ionization; state mapping and unresolved final products are
explicit in docs/ATOMIC_NETWORK.md. All three exact personal files are now
preserved outside Git in
`/Users/max/Library/Application Support/FRC/personal-atomic-inputs/` as well as
their original temporary paths. No provider tables/rates are bundled here.
