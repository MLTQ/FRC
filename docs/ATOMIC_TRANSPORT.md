# Spatial transport of the reacting moment state

The actual heavy-level and electron state now has conservative axial fluxes
for number, all three momentum components and total translational energy.
`AtomicTransport` composes those fluxes with the existing reacting/colliding
source. Densities use one fixed reference number density, not a cell-by-cell
normalization. See `plasma_bench/atomic_transport.md` and
`plasma_bench/species_euler_flux.md` for the equations and fixed-frame pressure
work. Binding energy travels with populations; produced photon energy is a
stationary local production ledger, not an advected radiation field.

The spatial method is minmod reconstruction, HLLC with HLL fallback and SSP-RK2.
Its numerical diffusion is explicit. No state, energy, charge or current is
repaired and there is no density floor. This component has periodic axial
boundaries. It does not yet provide electromagnetic feedback, radial geometry,
physical open ends, neutral collisions or a tracked vacuum interface.

## Complete independent transport and decay histories

`outputs/atomic_transport` contains all six completed histories, every accepted
state, analytic final values and integral budgets. Run 84798 completed in
132.214 s; source
`c50e5c4e5369839b839effff7cadf2a2b83175d7499478e1e41660b316383b58`.
The scientific PNG was inspected and the PDF is saved. Raw personal atomic
tables are not bundled. These are mathematical component fixtures, not lab
preparation or formation runs.

An initially charge-neutral electron/ion contact with constant separate pressures
translates through one full 1 m period in 10 us. Physical electron and argon
ion masses are retained. The exact charge and current are zero. A second case
translates an uncharged excited argon population through 1 cm in 40 ns while
it decays through a supplied single radiative channel. Its reference includes
the independently integrated local photon-production history.

| Case | Cells | Accepted steps | Final mean population L1 error |
|---|---:|---:|---:|
| Plasma contact | 32 | 393 | 0.011953813 |
| Plasma contact | 64 | 789 | 0.003950076 |
| Plasma contact | 128 | 1580 | 0.001146706 |
| Neutral decay | 16 | 54 | 0.031206598 |
| Neutral decay | 32 | 107 | 0.012049943 |
| Neutral decay | 64 | 214 | 0.003998365 |

Fixed acoustic CFL refines space and time together. Contact L1 errors decrease
by factors 3.03 and 3.44; limiting near smooth extrema prevents an assertion of
uniform second-order accuracy. Neutral photon maximum errors decrease from
0.0550731 to 0.0152174 to 0.00410183 eV per reference nucleus. All neutral
electron numbers and charges remain exactly zero. Contact local charge error
is at most 7.55e-15 reference charges. The contact integral energy drift is at
most 1.43e-14 eV m per reference density, with nucleus/mass errors below 1e-15 m.
Neutral energy drift is at most 2.23e-16 in the same energy-integral units.
These are observed integral budgets, not imposed conservation corrections.
Direct laboratory current reconstructed from every accepted extensive state
is at most 1.23e-9 A/m²; `current_readout.json` retains the per-grid values.

## Pressure work and the vacuum-front limitation

`outputs/atomic_neutral_expansion` contains three complete neutral slab
expansions through 9.850823 us (10.055 s wall time), source
`69e5b29ef71c26cb9bc2ee5a6f9374fe68dd35d243a47d5776b4c031d2058e46`.
The terminal summary reports success; the earlier process handle is expired.
All histories and exact rarefaction cell averages are saved; the scientific PNG
was inspected and PDF saved. The independent gamma=5/3 vacuum rarefaction is
specified in `examples/atomic_expansion.md`. Unlike a translating contact,
this exercises pressure work and accelerating flow.

| Cells | Steps | Density mean L1 error | Fraction of nuclei beyond exact front |
|---|---:|---:|---:|
| 64 | 25 | 0.005827129 | 0.000368898 |
| 128 | 72 | 0.004581425 | 0.000147640 |
| 256 | 208 | 0.002779996 | 0.000074776 |

The finest scaled-momentum/energy mean errors are 0.00108297 / 0.000583431;
all charged species stay exactly empty. Maximum nucleus and energy integral
drifts are 8.33e-17 m and 1.39e-17 eV m per reference density.

**This Eulerian flux does not preserve the exact free boundary.** Numerical
transport puts nonzero density in every final grid cell. On the finest grid,
the maximum density in cells wholly beyond the analytic front is 0.000967948
of the initial density. The leak decreases with refinement but must not be
renamed vacuum or clipped away. It disqualifies this scheme as a replacement
for the formation solver's tracked plasma/vacuum interface. No formation run
has adopted it; the existing exact exterior vacuum remains intact.

The next physical coupling must retain actual species fluxes and their charge,
momentum and energy in the electromagnetic/current-inertia equations. A local
population update that resets the moving solver's fixed carrier weights or
two entropies would discard those equations. Boundary tracking and physical
escape remain required alongside that coupling; these fixtures do not close
the full FRC goal or the atomic/moving integration issue.
