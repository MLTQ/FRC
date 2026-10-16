# Radiative recombination with electron, atomic and photon energy

The atomic source now evolves adjacent recombination alongside direct impact
ionization. This remains a fixed-volume component; it is not yet connected to
the moving pulse. `plasma_bench/radiative_recombination.md` states the equations,
units and assumptions; `examples/radiative_recombination.md` describes the
complete trajectories. Data provenance is in `data/atomic/recombination/README.md`.

For each recombination event the source removes the captured electron's mean
kinetic energy from electron heat, lowers the atomic binding inventory by the
successive ionization energy, and puts both contributions into emitted photons.
The mean capture energy is derived from the Maxwellian rate moment,
T[3/2 + d ln(alpha)/d ln(T)]. It does not assume captured electrons have the
bulk-average energy 3T/2. Thus remaining electrons may warm even as electron
thermal energy decreases. No state floors, clipping, energy restoration,
renormalization, unit tests or verification framework are introduced.

This requires immediate cascade to ground and escape of all photons. Neither
is established in the experiment. The present source excludes excited-state
storage, opacity/reabsorption, line cooling, DR, three-body recombination,
charge exchange, neutral motion and species momentum/current exchange. Fit
derivative accuracy is an additional limitation beyond rate-value accuracy.

## Complete runs (2026-10-16)

Run 3052 completed with exit 0 in 11.218 s. The frozen source hash is
`5523f224514cdfd8c3d182612e3ebd500f2d4a038f035f9c9fe80bac0098778a`;
RR input hash is
`014e2afc48aca561002ba642c6f41ddb28ad06ac8feb214b21d229af353956f6`.
`outputs/radiative_recombination` contains 15 complete RR/reference trajectories,
four matched direct-ionization-only controls, input/source copies, summary,
individual original-rate comparisons, and inspected PNG plus PDF.

The independent power-law recombination solution has maximum ion-fraction
errors 2.90594e-8, 1.84431e-9 and 1.15861e-10 for 50/100/200 steps over 30 ms.
Successive improvements are 15.76 and 15.92, consistent with RK4. The finest
maximum temperature error is 7.13e-11 eV and photon-energy error 1.904e-9
eV/nucleus. Thermal + atomic + photon energy residual is <=7.22e-15 eV/nucleus.
The analytic fixture specifies a power-law cross section; it is not an
experimental hydrogen recombination dataset.

All argon cases use fixed N=8.961429541213399e19 m^-3. This matches the initial
central pulse nucleus density, not the tube-averaged neutral fill. Finest runs:

| Initial state / prescribed heating | Duration | Final mean charge | Final Te, eV | Emitted RR energy, eV/nucleus |
| --- | --- | --- | --- | --- |
| All Ar+, 10 eV, no heat | 3.5 µs | 1.217345591 | 4.914237577 | 0.001592019 |
| 10% Ar+, 10 eV, no heat | 3.5 µs | 0.154005350 | 2.775031689 | 0.000019108 |
| 10% Ar+, 10 eV, 20 eV/nucleus over 1 µs | 3.5 µs | 0.974093388 | 4.111391792 | 0.000816021 |
| All Ar+, 1 eV, no heat | 1 ms | 0.970377511 | 1.013045356 | 0.492284246 |

For the first case, excluding RR gives final mean charge 1.217395187 and
Te=4.914273882 eV. The RR photon loss is small compared with about 6 eV/nucleus
stored in further ionization. This only compares those two included channels;
it says nothing about the size of missing excitation/line-radiation losses.

In the 1 ms case electron thermal energy falls from 1.5 to 1.474554646
eV/nucleus despite Te rising. This is selective capture of slow electrons,
not heating from binding release. It is not a physical afterglow prediction.

The final timestep halving in the three short cases changes populations by
<=3.10e-9, Te by <=1.40e-7 eV and photon energy by <=7.04e-12 eV/nucleus.
The cold case is already at floating-point sensitivity at all three levels.
Minimum population is zero; maximum nucleus-number residual is 2.89e-15.
Energy residual in unheated runs is <=5.29e-14 eV/nucleus; the finest prescribed
heating case gives 1.004e-8 eV/nucleus against the independent analytic integral.
No claim of exact convergence order is made for the heating pulse endpoint.

## Data limitations and next physical work

The author's original tabulated RR calculations supply 228 temperature points
over 12 charge stages. The maximum fit-versus-rounded-original difference is
3.802%, at Ar15+; stage maxima range from 0.744% to 3.802%. Four modern stages
(8,16,17,18) lack original rows in the downloaded master. This is a comparison
to an independent atomic calculation, not a laboratory validation, and does
not validate the energy moment obtained by differentiating fits.

More seriously for these examples, Ar+ and Ar2+ require older published
power-law fits whose validity intervals are unspecified. They supply about
99.25%, 99.9994%, 98.25% and 100% of time-integrated RR events in the four cases.
Their uncertainty is exposed, not eliminated by numerical convergence.
`event_content.json` integrates source event rates at saved states with the
trapezoidal rule. The cold run has only 2.923e-9 out-of-range impact-ionization
events per nucleus, 0.0366% of its integrated ionization events, despite a
100% *instantaneous* out-of-range fraction initially when all target ions are
Ar+ and Te is below that fit's 2 eV lower bound. Both facts matter.

Next: obtain defensible low-charge excited/metastable, DR and three-body data;
add excitation/line-energy storage and radiation with explicit density/opacity
limits; derive mass, momentum, current and energy transfers before moving
adoption. The four unchanged long pulse jobs continue their time/spatial
comparisons independently. The full FRC objective remains open.
