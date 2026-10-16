# Electron-impact excitation and resolved energy

Finite atomic level source with Maxwellian excitation/de-excitation, spontaneous
radiation and a finite electron thermal reservoir. It is independent of the
still-running moving pulse. Charge changes, neutrals' momentum, ion heat,
non-Maxwellian/finite-drift excitation, photon absorption and spatial transport
remain absent; the level source alone is not a plasma cooling closure.

For supplied effective collision strength Upsilon(T), excitation energy DeltaE,
and lower/upper weights g_l/g_u, the SI rate coefficients are

    q_up = C Upsilon/(g_l sqrt(T_K)) exp(-DeltaE/(kB T)),
    q_down = C Upsilon/(g_u sqrt(T_K)),
    C = sqrt(8/(pi*m_e))*pi*a0^2*Rydberg/sqrt(kB).

This is the Maxwellian cross-section integral with
sigma(E)=pi*a0^2*Rydberg*Omega(E)/(g_l E). C is approximately8.629e-12 in SI.
Detailed balance follows directly. Upsilon is interpolated linearly in log(T)
with nonnegative endpoint weights; source temperatures bound the model and no
extrapolation is allowed. The literal source table is neither fitted nor floored.

Population transfers use ne*q*f with both directions, retaining all supplied
pairs. Electron energy per modeled atom receives -DeltaE*(R_up-R_down).
Each spontaneous transition adds its level-energy difference to a separate
produced-photon energy reservoir. Atomic excitation energy is sum(E_i*f_i).
The population/electron/photon generator therefore has exact particle and total
energy left invariants. Electron energy is evolved as its own state variable;
it is never reconstructed from a conservation residual or repaired afterward.

The neutral input also contains a ground-state ionization S row. A fixed-charge
excitation component must explicitly pass `exclude_processes=("S",)`; otherwise
construction rejects. This is a declared component limit for comparison, not a
replacement of the existing charge-changing source. Integrating charge changes
requires one consistent set of rates/thresholds, electron-number and energy
updates, without counting ionization twice. Ar II's inspected input has no
additional process rows. Missing/unmodeled processes are never certified zero.

`step` uses an exponential half-step at initial electron temperature to predict
the midpoint temperature, then a full exponential at that temperature. This is
second-order explicit exponential midpoint for the nonlinear reservoir feedback,
not an exact solution of the nonlinear source. No floor or positivity repair is
used; invalid states or out-of-range temperatures reject. Complete time refinement
and independent finite-reservoir Boltzmann equilibrium are required measurements.

`radiative_path` optionally replaces energies and radiative A values with an
explicit independently supplied record, only when ordered degeneracies match
and every energy differs by <=1 cm^-1. The caller must disclose that replacement.
The Ar I application uses current NIST data matching the existing15 levels;
source A values differ substantially and are not silently treated as current
NIST values. Ar II source states are LS terms: source xj is a statistical-weight
encoding, not actual fine-structure J. Its term-averaged radiation is not a
resolved line profile. The constructor currently assumes argon's mass39.948u.

Raw OPEN-ADAS tables are explicit personal local inputs, kept outside the repo
and omitted from source/output archives; this module does not download or bundle
them. Input paths and hashes identify which physical rates were used. Provider
restrictions on bundling/publication are not removed by this parser.
