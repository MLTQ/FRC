# Joint atomic populations, charge and electron energy

Combines existing full ground charge ladders with supplied resolved excitation
networks. Every nucleus has one population coordinate: q=0 and q=1 may have
many levels, while higher argon charges currently have ground coordinates only.
This removes the separate fixed electron bath from the excitation component.
It is a homogeneous source, not yet the moving finite-mass/current simulation.

State is `[populations, electrons/nucleus, electron energy, produced photons,
integrated heat]`; energies are eV/nucleus. Electron number and energy
are independently integrated. Quasineutrality is a conserved source invariant,
not restored by overwriting electron number with mean charge. Temperature is
Ee/(1.5*electron_number); ne=N*electron_number. Exact zero carriers produce no
impact events or thermal heating; there is no artificial seed or floor.
Spontaneous decay still evolves an excited neutral population with exactly zero
electrons. It is assembled independently of temperature-dependent collisions.
The first source archive predates this zero-carrier correction; its completed
positive-electron histories remain valid records of that older source.

Atomic energy epsilon includes cumulative independent NIST charge thresholds
plus supplied excitation energies. For an ionization event a->b the population
transfers once, electron number increases once, and electron energy loses
epsilon_b-epsilon_a. Excitation/de-excitation keeps charge fixed and transfers
the level spacing between atom and electrons. RR removes an electron, removes
the fitted Maxwellian mean capture energy from electron heat and produces
binding-release plus capture energy as photons. All channels are added to one
generator; nucleus count, electron_number-sum(q*f), and
Ee+sum(epsilon*f)+photons-integrated_heat are exact generator invariants.

The matrix is evaluated at a predicted exponential midpoint, updating both ne
and Te. All populations/number/energies evolve together. Invalid stages or
endpoints reject; no clipping, normalization, charge fix or energy repair.
Sparse exponential action preserves the disconnected zero-population blocks.
The earlier dense Padé exponential generated negative populations around1e-18
in an entirely unpopulated Ar II block during zero-electron neutral decay; that
rejected state and accepted prefix are preserved, not rounded away. An appended
constant coordinate represents the prescribed affine heat source only during
the exponential action; it is not part of the evolved physical state.
The numerical time method is second order, not an exact nonlinear solution.
Heat is a nonnegative prescribed source integrated in its own coordinate; an
independent analytic heat integral measures its time-discretization error.

## Explicit data/model limitations

Resolved ionization mappings are (charge, local level index, ADF07 row). They
replace, rather than duplicate, a ground fit when the origin matches. The
argon example identifies 2P, 4D and 4F LS terms explicitly; source metastable
IDs must never be used directly as array indices. Source effective potentials
remain separate from the model's physical threshold, with discrepancy metrics.

Final ADF07 code0 is unresolved. This preliminary joint source collapses such
products to the next ground coordinate. Final excited branching/storage and
its effect on electron energy are unavailable, not established zero. RR likewise
retains the older immediate ground-cascade approximation; it does not populate
the resolved bound states. Only supplied initial ground RR fits are used;
excited parent RR is omitted explicitly. These assumptions preclude claiming a
complete collisional-radiative model or net physical cooling.

Unsupplied excited-state ionization is omitted and its occupied population is
reported. It is not assigned a ground-state rate, Boltzmann population, fitted
threshold correction or a demonstrated physical zero. Missing higher excitation,
dielectronic/three-body recombination, neutrals/momentum/current transfers,
finite drift, photon absorption/transport and physical walls remain required.
Produced photons are not established escape. Existing data range and legacy
RR uncertainties remain visible. Do not install this as an unqualified plasma
energy sink. Complete joint trajectories quantify the coupled consequences and
the unsupported state populations before further physical adoption.

The generic solver embeds no OPEN-ADAS tables; all personal restricted inputs
are supplied outside the repository. Ground-only analytic DI/RR limits use the
existing independently derived references, without a new test framework.
