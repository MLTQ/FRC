# Explicit ground-argon cross-section replacement

`ArgonGroundCrossSection` supplies the existing `AtomicNetwork` coefficient
interface, replacing only Ar(ground) -> Ar+ with the integral of a private
LXCat `IONIZATION / Ar <-> Ar^+` record. All other charge rates, physical NIST
atomic energies and charge coordinates stay identical to the base provider.
The base fit and replacement never contribute simultaneously to this reaction.

The source's threshold is retained for discrepancy diagnostics; the atomic
electron-energy loss continues to use the independent physical ionization
energy, exactly as in the existing fit-based model. No source-energy-axis shift,
threshold fitting, fitted scale, state floor or energy/charge repair is made.
The `<->` source label does not invent an inverse recombination process.

The finite energy range is integrated without extrapolation. Since the table
does not prescribe a validated Maxwellian-temperature interval, the returned
range mask is false for the replaced channel. Its events therefore contribute
to the existing outside-fit-range diagnostic; this means unspecified range,
not that all sampled collision energies lie outside the table. The example
separately reports the electron flux fraction above its upper energy limit.

This input does not supply excited-parent ionization or final excited-product
branching. Existing ground-product and radiative-recombination assumptions stay
explicit. The adapter does not alter moving pulse inputs or adopt an atomic
source into the formation solver. Source tables stay outside Git.
