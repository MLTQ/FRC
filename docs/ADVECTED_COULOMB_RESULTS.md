# Moving nonuniform plasma with electron-ion heat exchange

Three full nonlinear runs completed one 100 us transit, with all fluid,
species, current and field equations active. The plasma moves at 2400 m/s
through a periodic 0.24 m cylinder of radius 0.035 m. Initial n=1e20 m^-3,
Ar+, both gamma=5/3, Bz=5 mT, Te=5+0.3 cos(2 pi z/L) eV and Ti=6-Te.
No resistivity, numerical damping, filtering or temperature floor is used.

The independent reduced-model solution is local Coulomb equilibration along
z-vt characteristics, with constant n, total pressure, velocity and B, and J=0.
It uses a separately integrated scalar temperature clock, not the production
RHS. This benchmark omits conduction and is not a complete laboratory model.
Full equations, boundaries and numerical approximations are documented in
examples/advected_coulomb_contact.md and the Coulomb model companion.

## Measured endpoints

The radial spaces remain F4/B8 because the exact state is radially uniform.
The table uses independent order-20 physical reconstruction. Temperature
errors are volume RMS; energy changes reference the actual initial state.

| Axial points / steps | Te error [eV] | Ti error [eV] | Relative energy change | Relative n RMS error |
|---|---:|---:|---:|---:|
| 16 / 80 | 1.70448e-5 | 4.93789e-5 | 6.28978e-10 | 8.65613e-6 |
| 32 / 80 | 3.41121e-7 | 8.85853e-7 | 6.29050e-10 | 1.52975e-7 |
| 32 / 160 | 3.25165e-7 | 8.67153e-7 | 3.93322e-11 | 1.51974e-7 |

Doubling axial resolution at fixed time step reduces Te/Ti error by 50.0/55.7
times. Halving the time step on the fine mesh changes those small temperature
errors little, but reduces energy drift by 15.99 times and entropy-change error
from 1.97884e-9 to 1.23659e-10 kB per initial particle. This separates the
remaining spatial temperature error from the dominant temporal energy error.
Maximum particle-number change over all runs is 3.12e-17 relative.

Small unwanted currents and fields do not improve with time refinement:
the 32/160 run reaches 7.85273e-3 A/m² and 1.08017e-11 T. Its maximum
relative current drift divided by the combined thermal speed is 6.48932e-10.
These departures are tiny in absolute magnitude but are not claimed converged.
This zero-current reference does not close the magnetic/current problem
frc-7mm.7 or justify the rejected hyperdiffusion candidate.

## Evidence and scope

Raw directories are outputs/advected_coulomb_contact_z16_t80,
outputs/advected_coulomb_contact_z32_t80 and
outputs/advected_coulomb_contact_z32_t160. Each contains the source closure,
initial state, all accepted endpoint and stage arrays, and solver history.
Actual sessions 38573, 39308 and 91750 exited 0 and were collected. Elapsed
times were 86.11, 132.70 and 273.04 seconds; peak RSS 276873216,
416514048 and 403505152 bytes. No worker remains active for this comparison.
No unit tests or new verification framework were created or run.

This establishes a bounded moving collisional component. Exact vacuum,
physical coil/end boundaries, conduction, ionization and credible formation
remain open. The next simulation couples temperature-dependent resistance,
magnetic decay and Joule heat under explicit weak-field transport assumptions.
