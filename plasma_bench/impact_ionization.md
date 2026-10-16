# Electron-impact ionization with atomic energy storage

An opt-in atomic source component, not yet connected to the moving pulse.
`ImpactIonization` accepts a complete charge ladder with independent physical
ionization energies and the published Voronov direct-ionization rate fits.
All internal energies are joules and densities are m^-3; rates are m³/s.
`from_file` reads an explicit data record and retains its hash. Argon data and
primary-source attribution are in data/atomic/argon/README.md.

For fixed nucleus density N, charge fractions f_z and quasineutral electrons,
ne=N*sum(z*f_z). Temperature is derived from evolving electron thermal energy
per nucleus e_e=(3/2)*sum(z*f_z)*T_e, with T_e in joules. Event rates per nucleus
are r_z=ne*f_z*k_z(T_e). Adjacent charge states lose/gain exactly those events:

    df_z/dt = r_(z-1)-r_z
    de_e/dt = heating_per_nucleus - sum(chi_z*r_z).

The atomic energy per nucleus is sum(f_z*epsilon_z), with
epsilon_z=sum_(j<z) chi_j. Thus electron thermal plus stored atomic energy
changes only by declared heating. The larger electron population also cools
through the temperature definition; no extra thermal-energy sink is counted.
There is no temperature/abundance floor, renormalization or state repair.
Negative trial populations/energy are rejected; a suitable integrator must
resolve reaction times. No artificial seed is introduced: exactly neutral,
zero-electron-energy gas has zero electron-impact source, and thermal electron
heating cannot be defined there without a preparation/seed-electron model.

Rate fits assume a Maxwellian electron distribution and ground-state targets.
Every coefficient evaluation returns its published temperature-range mask.
Source diagnostics quantify event fractions outside that range and involving
fit thresholds differing by more than5% from physical ionization energies.
These are review indicators, not permission to extrapolate or validity proofs.
The raw coefficients are never modified to conceal either limitation.

This source contains direct ionization only. It has no recombination, excited
or metastable populations, radiative cooling, charge exchange, particle
transport or recoil/three-fluid momentum exchange. A homogeneous source-only
trajectory is an equation/component study, not a gas-preparation prediction.
Before moving adoption, all applicable processes and the new particle/mass/
current/energy transfers must be coupled explicitly with exact exterior vacuum.
