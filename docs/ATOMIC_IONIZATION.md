# Atomic charge evolution and electron energy

The moving pulse currently assumes a fixed singly-ionized inventory. This
omits both changing electron number and the thermal energy spent on further
ionization. New `plasma_bench/impact_ionization.py` evolves a complete direct
electron-impact charge ladder with independent atomic-energy bookkeeping.
It is a source component; moving-pulse coupling remains unfinished.

## Physical contract and data

For fixed nucleus density N and charge fractions f_z, ne=N*sum(z*f_z).
With electron thermal energy per nucleus e_e, Te=e_e/[(3/2)*sum(z*f_z)]
in joules. Events r_z=ne*f_z*k_z(Te) move one nucleus from z to z+1 and
create one electron. Each event removes the physical ionization energy chi_z
from electron thermal energy and adds it to atomic stored energy. Temperature
also falls as that thermal energy is shared among more electrons. These two
effects are not counted twice. The invariant is electron thermal plus atomic
energy, changing only by the prescribed external heating in this source study.

Rates use [Voronov's published fit table](https://www.pa.uky.edu/~verner/col.html)
for ground-state direct ionization, with its original temperature limits.
Ground-state energies and uncertainty/estimate annotations come separately
from [NIST ASD](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html).
All18 argon transitions, raw source files, URLs and hashes are saved in
data/atomic/argon. In particular the archival Ar3+ rate-fit threshold52.3eV
differs from NIST's59.58eV. Neither the fit nor the energy is silently changed
to agree with the other; reaction-weighted content of that discrepancy is
reported. The present data require further review before broad physical use.

The rate formulas assume Maxwellian electrons. No recombination, excitation,
metastables, radiative cooling, opacity, transport, recoil or new-species
momentum exchange is included. There is no artificial electron seed, population
floor, clipping, renormalization or energy repair. Negative intermediate
populations/energy reject the attempted calculation. Exactly neutral gas with
zero electron energy has zero electron-impact source; preparation needs a
physical seed/ionization mechanism. No exterior atmosphere was introduced.

## Complete source calculations, 2026-10-16

Run60848 exited0 in3.122s. All12 complete trajectories, source archive and
input-data copies are in outputs/impact_ionization. Source:
d418d918fa571267614effdc7962f24f6f397567f80527a709594b27e4099a85.
Argon record SHA256:
6f90efc5cc8267d63a5765d2d2e0f0529595c337710ee1e58117452651af93f0.
The scientific PNG was inspected and PDF saved. No unit tests were run.

An independent hydrogen reference eliminates thermal energy analytically and
integrates elapsed time as a scalar quadrature over ionized fraction. It does
not call the production source/rate implementation. Complete20us RK4 runs
with100/200/400steps have maximum ion-fraction errors4.46844e-6 /2.92686e-7 /
1.87101e-8 and temperature errors3.39457e-4 /2.35477e-5 /1.52840e-6eV.
This is approximately fourth-order improvement. The largest energy drift is
1.44249e-14eV per nucleus. It verifies the source evolution against an
independent mathematical solution, not the experimental accuracy of a rate fit.

Argon calculations use fixed N=8.961429541213399e19m^-3, the pulse's initial
central nucleus density, and initial Te=10eV. All19 charge states evolve for
3.5us with400/800/1600steps. This is a uniform parcel, not the tube's mean fill
or a material trajectory from the moving pulse. The heated case deposits
20eV per nucleus over1us using a prescribed sin² thermal-power envelope;
it does not define an electrical source or experimental preparation method.

| Initial parcel / source | Final mean charge | Final Te | Extra stored atomic energy | Final neutral fraction |
|---|---:|---:|---:|---:|
| All Ar+, no heating |1.217395|4.914274eV|6.026080eV/nucleus|0|
| 10% Ar+ seed, no heating |0.154006|2.775028eV|0.858942eV/nucleus|0.846653|
| 10% Ar+ seed, prescribed20eV heating |0.974123|4.111383eV|15.492513eV/nucleus|0.168205|

The first case ends with21.4417% Ar2+ and0.14886% Ar3+. Mean charge does not
equal the fraction of ionized nuclei; the heated case still has16.82% neutrals
despite mean charge near1. These are consequences of the declared source-only
equations, not final predictions for the lab or the pulse with omitted physics.

Across all argon histories the minimum population is0, nucleus-number error
<=2.78e-15. Halving the last timestep changes populations by at most3.10e-9
and Te by1.40e-7eV. The heated energy residual against an independent analytic
heating integral decreases1.54224e-6 ->1.44576e-7 ->1.00400e-8eV per nucleus;
its pulse endpoint does not coincide with every uniform grid. Unheated energy
residuals remain below2.89e-14eV per nucleus. These are numerical measurements,
not corrections applied to the state.

Maximum fraction of reaction events using an extrapolated rate is3.13e-14
in these cases; maximum event fraction involving the >5% threshold discrepancy
is8.09e-6. Those small fractions do not validate the Maxwellian distribution,
omitted metastables/radiation/recombination, or other densities/temperatures.

## Remaining physical work

frc-transport.5 stays in progress. Add applicable recombination, excitation and
radiative energy channels from defensible atomic data, then couple the evolving
charge/neutral populations to particle, finite-mass momentum, current and
species energy equations. Newly ionized material cannot silently inherit a
velocity or fixed-Z two-fluid representation. Keep nonlocal heat and physical
edge/end particle/energy budgets active underfrc-transport.4. The four existing
pulse refinements retain their unchanged sources and cannot be retroactively
treated as simulations containing these new atomic processes.
