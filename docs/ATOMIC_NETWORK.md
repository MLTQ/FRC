# Coupled atomic levels, charge and electron energy

The new source evolves resolved excitation and charge changes in the same
finite electron reservoir. Electron number is independently advanced along
with the populations, thermal energy and produced radiation. This replaces the
separate fixed-charge excitation bath for these homogeneous source studies.
Moving fluid/current coupling and complete radiation transport remain unfinished.

Final run 55867 exited 0 in 89.663 s: 12 full argon trajectories and seven
independent-limit trajectories. `outputs/atomic_network_sparse` contains every
accepted state, exact source archive, input identities, summary and inspected
PNG/PDF. Source is
`25eb4512063bdfb10f3395c52d8f85b31f04d6d976854dd341de8b7e26448a33`.
No unit tests or new verification framework were used.

## State, equations and physical assumptions

There are 72 nucleus-population coordinates: 15 Ar I levels, 40 Ar II LS terms,
and one ground coordinate for each charge 2 through 18. For nucleus density N,
the separately evolved electron number x gives ne=N*x and Te=Ee/(1.5*x), with
energies in eV/nucleus. Atomic energy epsilon contains the independent cumulative
NIST ionization energies plus each supplied excitation energy.

An ionization event a->b transfers one nucleus, creates one electron and removes
epsilon_b-epsilon_a from electron energy. Excitation/de-excitation exchanges the
actual level spacing with electron heat at fixed charge. RR removes an electron,
subtracts its fitted Maxwellian capture-energy moment from electron heat, and
produces that energy plus the atomic binding release as radiation. The source
has three simultaneous invariants: nucleus count; x-sum(q*f); and
Ee+sum(epsilon*f)+produced_photons-integrated_heating. No coordinate is overwritten
to enforce an invariant. All supplied forward/reverse excitation channels remain.

The new [ADF07 Ar II input](https://open.adas.ac.uk/element/detail/adf07/szd97%5D%5Bar/szd97%5D%5Bar_ar1.dat)
supplies three initial LS terms, 2P, 4D and 4F, over 0.6–4000 eV. Source
metastable IDs 1/2/3 map to excitation indices 0/2/5 by configuration,
multiplicity, orbital L and statistical weight. The ground Ar II table replaces
its old ground fit; it is not added a second time. Positive log-log interpolation
stays within the input interval. Its effective potentials are retained separately
from the model's energy differences. The source comments describe Griffin's
1996 calculation and partial normalization to experimental data; this does not
establish experimental accuracy of the combined plasma model.

Final metastable code 0 is unresolved. The provisional source collapses those
products to the next ground coordinate; final excited branching, storage and
its electron-energy cost are unavailable. Total RR fits likewise retain the
earlier immediate-ground-cascade assumption and do not populate resolved levels.
Only initial ground RR fits are supplied. Unsupplied excited ionization and
excited-parent RR are omitted explicitly, with occupied-population diagnostics;
no ground rate is assigned to them. Neutral ADF04's S row is excluded in favor
of the one existing ground DI fit, avoiding duplicate ionization.

Maxwellian electrons, fixed nuclei, prescribed homogeneous heating, missing
higher excitation, DR/three-body recombination, recoil and neutral/current
momentum transfer limit this source. Produced photons are not established
escape or net cooling. The earlier opacity and finite-frequency radiation
limitations remain. These assumptions must be resolved before a complete
moving-plasma energy closure is claimed.

## Time integration and a preserved numerical failure

A half exponential predicts the midpoint ne/Te; a full exponential at that
midpoint advances all physical coordinates together. The method is second
order. Prescribed heat uses an affine source coordinate only during the
exponential calculation, plus an independently evolved heat integral.
Negative/nonfinite populations or energy reject without clipping or repair.

Initial run 75960 completed 18 trajectories in 19.195 s, source
`6102b74ecf5a2acd991dae9b88ef6676e88c8bf9301a18fe89c47e9a16185ce0`,
in `outputs/atomic_network`. Inspection then exposed an incorrect zero-electron
early return that also stopped neutral spontaneous decay. That path was corrected.

Intermediate run 52694 exited 1 in `outputs/atomic_network_zero_carrier` after
the six analytic charge trajectories. Its first neutral midpoint contained
negative populations down to -1.41828e-18 in the disconnected, initially empty
Ar II block. Dense Padé factorization caused cross-block roundoff; the constant
coordinate was exactly 1. The rejected array, initial state and terminal
diagnosis from run 19428 are preserved in that output directory. The actual
generator contains provider radiative rates and is retained only in the private
personal-input directory, with its path/hash in `neutral_failure.json`. Source:
`76c88a8a693a016f5e246a52245f4f9cfd32289732e9c6d34bb4a28ea148b347`.

The final method uses sparse exponential action, retaining the disconnected
zero populations. It does not round away the negative state from the earlier
attempt. Spontaneous radiation evolves independently of electron collisions.

## Complete argon histories

All use N=8.961429541213399e19 m^-3, initial Te=10 eV and duration 3.5 us.
The seeded cases begin 10% singly ionized and 90% neutral. The heated seed
receives 20 eV/nucleus through a 1 us sin-squared pulse. Joint runs use
700/1400 steps (5/2.5 ns); controls use 1400 steps.

| Case | Model | Final mean charge | Final Te, eV |
|---|---|---:|---:|
| Initially all Ar+ | Joint |1.04927679|2.31777975|
| Initially all Ar+ | Joint without two excited ionization channels |1.03185547|2.31209116|
| Initially all Ar+ | Older ground-only |1.21734553|4.91423878|
| 10% seed | Joint |0.13564697|2.38680445|
| 10% seed | Joint without two excited ionization channels |0.13560006|2.38714661|
| 10% seed | Older ground-only |0.15400533|2.77503358|
| Heated 10% seed | Joint |0.68390234|2.53628064|
| Heated 10% seed | Joint without two excited ionization channels |0.67118062|2.52826197|
| Heated 10% seed | Older ground-only |0.97409333|4.11139122|

The old ground-only control also uses a different ground Ar II rate source;
its difference must not be attributed solely to excitation. The excited-channel
ablation holds that ground source fixed. Its two new channels account for
36.5718%, 0.1923% and 2.5702% of integrated ionization events in the three
joint cases. These event integrals use saved-time trapezoidal quadrature.

Maximum occupied fractions without supplied excited ionization are respectively
3.1480%, 0.3327% and 2.8287% of nuclei. Missing excited-parent RR populations
reach 5.9964%, 0.1372% and 3.5590%. The ground-only control's zero on that plot
means it has no excited coordinates, not that all physical channels are known.
These are populations, not bounds on missing reaction power. New bead
frc-transport.5.4 tracks these data/branching requirements.

Across the three joint time-halving comparisons, maximum population L1 change
is 1.72801e-5, charge change 8.79874e-6 and temperature change 0.00104980 eV.
Across all 12 argon runs, nucleus-count error is <=2.72005e-13, charge-budget
error <=1.40554e-13 and energy-invariant error <=4.89564e-12 eV/nucleus. The
separate analytic-heating residual is larger: 1.30903e-4 eV coarse and
3.27252e-5 eV fine in the heated case, a fourfold decrease. Keeping these two
budgets distinct prevents a conserved numerically integrated heat source from
being mistaken for exact forcing accuracy. Legacy low-charge RR range is still
unknown; its event fraction reaches 1.0. The 7.28 eV archival high-charge fit
threshold discrepancy is retained, with reaction-weighted content in the summary.

## Independent mathematical limits

The existing scalar hydrogen ionization quadrature gives maximum charge errors
3.10097e-4, 7.33802e-5 and 1.78333e-5 at 100/200/400 steps through 20 us.
The existing exact finite-reservoir power-law recombination solution gives
8.96384e-5, 2.18238e-5 and 5.38425e-6 at 50/100/200 steps through 30 ms.
Both show the new method's expected second-order trend. These are independent
physical limits, not full atomic experimental validation.

The 40 ns neutral resonance decay starts with exactly zero electrons and
electron energy. Upper population matches exp(-A*t) within 4.44089e-16 and
photon energy matches DeltaE*(1-exp(-A*t)) within 7.99361e-14 eV. Electron number
and thermal energy remain exactly zero without an artificial seed.

## Personal input preservation and next work

Exact restricted inputs are also preserved outside the repository at
`/Users/max/Library/Application Support/FRC/personal-atomic-inputs/` so temporary
cleanup will not erase them. Original /private/tmp paths remain in the immutable
run records. The new `szd97_ar_ar1.dat` SHA256 is
`6129b96ad78059983131baa7a9f228ae5365b06dd03c3d5064f1ed29c096115f`;
the two excitation hashes are unchanged from COLLISIONAL_LEVELS.md. No raw
provider data or reconstructed rates are bundled with code/output archives.
The [provider terms](https://open.adas.ac.uk/terms-and-conditions) still govern
future redistribution. No pending LXCat terms were accepted.

Next resolve populated missing channels and final-state energies/branching,
then connect atom/photon absorption and moving particle/momentum/current/energy
sources. Keep frc-transport.5.3/.5.4, parent .5 and all original FRC gates open.
