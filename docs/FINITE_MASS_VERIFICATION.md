# Moving finite-mass verification and remaining pressure error

The full suite passes **328 tests in 120.97 s**. All **21 conservative-rate
cases** and **18 full coupled trajectories** finish and pass archive/configuration,
physical accounting and independent-reference audits. This is a positive-density
bounded reference, **not an accurate FRC formation solver**. Global conservation
passes while local species-pressure accuracy remains inadequate.

The equations, exact masses, numerical diffusion, state partition, unsplit RK4,
boundaries and failures are specified in [FINITE_MASS_NUMERICS](FINITE_MASS_NUMERICS.md).
The scientific figure was visually inspected:

![Measured rate, trajectory and reconstruction errors](../outputs/finite_mass_evolution/accuracy.png)

## Independent references and protocol

The nonlinear manufactured reference has varying density, all velocity/field
components, two pressures and forces/heat evaluated directly from both parent
species. Initial conservative quantities and exact time derivatives use separate
physical r,z quadrature. Initial total-energy representation error is measured;
no accepted state is reinitialized to the reference. Parameters include physical
Ar+ and electron masses, R=.035 m, L=.24 m, B0=.005 T, n0=1e18 or 1e20 m^-3,
gamma_e=5/3, gamma_i=1.4 and a 1e7 s^-1 reference frequency.

The unforced wave oracle uses an independently eliminated finite-mass dispersion
relation and Bessel/Fourier spatial derivatives. All three branches satisfy the
original linear species/Maxwell equations; physical argon, me=0 and equal-mass
limits are tested. The n0=1e18 branch frequencies are approximately 2807.94,
545649 and 53924707 s^-1. The slow/medium branches have **rate** evidence only;
their long periods have not been time-integrated by this short trajectory study.

The 21 rates comprise five manufactured grids nr8→128, four original-density
grids nr8→64, three branches at nr8/16/32, slow branch nr64 and two tenfold weaker
nr32 modes. All meshes have nz=2 nr. The 18 trajectories comprise manufactured,
original-density and unforced nonlinear runs at nr8/16/32; three manufactured
time controls; fast waves at nr8/16/32; two fast time controls; and a weaker
nr32 wave. Manufactured/unforced duration is 50 ns, from 23 to 73 ns. Fast waves
run 29.1294366 ns, one quarter period. Base maximum dt is 1 ns, with explicit
2/.5/.25 ns controls as applicable. Exact steps and final fractional steps remain
in every history. Positive accepted states do not establish long-time stability.

## Complete-rate accuracy

Errors below are relative volume L2, with dual-volume weights for nodal pe_t.
Conservative and species rates use independent analytic targets; pointwise
maxima and first-annulus norms are retained in the JSON records.

| nr | Manufactured mass rate | Manufactured pe rate | Manufactured pi rate | Slow-wave pi rate |
|---:|---:|---:|---:|---:|
| 8 | .0704111 | .577037 | 5.34159 | 921.429 |
| 16 | .0232507 | .191818 | 2.00920 | 330.783 |
| 32 | .00812106 | .0649571 | .643218 | 117.248 |
| 64 | .00273643 | .0223374 | .233487 | 41.4570 |
| 128 | .000861842 | .00774879 | .143364 | not run |

At nr128 the manufactured radial/axial/angular momentum-rate errors are
.00143595/.000367828/.000416125; total-energy-rate error is .000529452.
The pi-rate **maximum does not decrease** from nr64 to nr128: about 30535.9 to
31209.6 Pa/s, despite decreasing volume error. The finest first-annulus pi-rate
error is .0125060. Pointwise consistency remains a separate unresolved gate.
Original-density pi-rate error reaches .0963559 at nr64. Medium/fast nr32 pi-rate
errors are .0283229/.454781. Tenfold smaller amplitude leaves slow/fast errors
essentially unchanged, so these errors are not explained by nonlinear amplitude.

## Full evolution and small-signal accuracy

Manufactured ion-pressure errors retain two denominators: full final pressure
and the true pressure change. Initial pressure representation error is near
roundoff, while initial total energy has a measured spatial representation error.

| nr | n18: full pressure error | n18: error in physical change | n20: full pressure error | n20: error in physical change |
|---:|---:|---:|---:|---:|
| 8 | .0501003 | 10.9270 | .0368497 | 8.03701 |
| 16 | .0138736 | 3.00257 | .00964106 | 2.08655 |
| 32 | .00362760 | .783580 | .00250372 | .540816 |

At n18/nr32, pi RMS/max/axis RMS errors are .000613613/.00203580/.000421504 Pa.
Electron-pressure error is .000269127 of full pressure but .0585963 of its
physical change. Conservative mass and total-energy change errors are
.00701935 and .00158038. The full pressure background must not hide the poor
resolution of its small physical evolution.

Fast-wave errors are normalized to the physical **oscillation**, not B0 or the
equilibrium pressures. Independent physical B/J and species pressure quadrature
includes their declared spatial representations.

| nr | Magnetic field | Electron pressure | Ion pressure |
|---:|---:|---:|---:|
| 8 | .158332 | .140771 | 15.4730 |
| 16 | .0795212 | .0384884 | 3.24947 |
| 32 | .0398052 | .0105524 | .595516 |

The finest current/density/velocity errors are .0203756/.0452679/.0400012.
Finest pi RMS error is 1.69434e-9 Pa. A tenfold weaker wave gives pi error
.595524, preserving the spatial limitation. Unforced nonlinear trajectories
have conservation and positive-state evidence; they have no analytic final-state
solution and are not labeled independently accurate physical trajectories.

## Time, conservation and entropy

Manufactured magnetic successive-halving ratios are 15.4502 and 15.3697;
fast-wave magnetic/current ratios are 15.4801/15.4703. Species pressure does
**not** show clean fourth-order convergence: manufactured pe ratios are
5.19414/3.30645 and pi 4.15122/3.14899; fast pe/pi ratios are 4.50654/3.18959.
Limiter switching, spatial coupling and small-signal recovery need further study;
RK4's formal order is not a measured universal order. Yet the changes from base
to half dt are far smaller than the spatial error: manufactured nr16 pi differs
by 1.08224e-7 Pa RMS versus 2.34673e-3 Pa error; fast nr16 differs by 4.09506e-12
Pa versus 9.24524e-9 Pa error. Smaller dt alone does not resolve this failure.

Maximum accepted-history relative mass drift is 4.44089e-16; total-energy budget
residual is 3.25241e-16. Independent physical component-energy replay differs
from stored total energy by at most 2.55e-14 relative. Manufactured parent
total-energy change independently checks source work; maximum absolute source
quadrature discrepancy is 5.07845e-12 J. Unforced axial/angular momenta and exact
shared-flux budgets are checked independently of near-zero net-wave sums.
Minimum accepted pe/pi are .734949/.138485 Pa. Full-E rate Faraday replay is
at most 3.48462e-11; maximum trajectory linear residual is 1.54128e-11 at the
unchanged 1e-10 tolerance. None of these algebraic properties proves accuracy.

Unforced physical electron log entropy decreases. With
D_s=(gamma_s-1) delta S_s/(N kB), electron values at nr8/16/32 are
-5.13023e-6/-1.49542e-6/-3.95357e-7; ion values are
.00404402/.000486844/.0000691347. Errors decrease but there is no entropy
inequality or verified shock heating partition. Physical r,z entropy quadrature
is distinct from nodal diagnostics. Forced entropy includes explicit source heat.

## Slow-wave cause and tested candidate

At n18/nr32, mass-rate error is 28.7004 times the physical rate. Its central part
contributes error 28.6816; the dissipative part has norm 1.03981. At nr64 these
are 10.1267/10.1204/.354951. Their scalar norms do not add. n20 gives the same
central-reconstruction dominance. The near-cancellation of large radial/axial
compression is lost by the limited second-order normal momentum traces.

Independent cubic reconstruction of true annular momentum averages improves
the **linear** continuity error to 1.93297/.197404/.0202495/.00194501/.000179355
at nr8/16/32/64/128. Polynomial exactness through axis/wall and smooth periodic
convergence have independent tests. Actual finite-amplitude initial momentum
instead reaches .00311608, exposing its nonlinear truncation plateau against
the linear oracle. This is a diagnostic, not an adopted solver or performance
claim. Consistent full momentum/enthalpy/kinetic transport remains to implement.

## Evidence and reproducibility

Raw initial rates, full fields and failed-accuracy diagnosis are preserved in
`outputs/finite_mass_initial_attempt`, source
`8937eefb61c197962881e7d844bee6d9f5dfce45aff93025c1c35e00dba63ea7`.
The 18 trajectories, all accepted histories, physical audits, complete final
fields, time comparisons, diagnostic arrays and figure are in
`outputs/finite_mass_evolution/report.json`, source
`3edae52928f5cdaa20ca240f11d6fa729a137be327c9b1ffc64f09c6440f42fb`.
Only the provenance driver's added evolution roots distinguish the old source
graph; all baseline production sources match. The report checks both identities,
every source/archive/record/NPZ hash and exact configuration. Helper/test/source
snapshots preserve the report's independent calculations. Prior moving-current
d3984fb... and current/thermal 80fd45b... fingerprints are unchanged.

Order-8 versus order-12 norm changes are at most 2.87874e-10 absolute in relative
L2 units. The manufactured axis heat evaluation at r/R=1e-10 agrees with
1e-12 within 5.63e-15 of its maximum scale; this reference limit is not a floor
in production physics. Rejected-step tests retain exact state, time and budgets.

Run the read-only study report with `python -m examples.finite_mass_report`
inside the project environment. Exact historical rate regeneration requires
its preserved source snapshot; do not overwrite or relabel it using a changed
driver. Fresh trajectory cases can be selected through
`python -m examples.finite_mass_evolution --group <name> --directory <new-path>`.
The complete test log and final checkpoint hashes are retained with the study.

Next improve shared spatial transport and repeat pressure/axis/maximum/time
checks, then physical collision/heating/transport and exact vacuum integration.
The full FRC goal remains open; no equilibrium, confinement or 3-D stability
claim follows from these bounded fixtures.
