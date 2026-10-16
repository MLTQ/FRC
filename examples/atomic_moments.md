# Full reactive-moment source histories

Runs the existing resolved argon network with separately evolved heavy-level
and electron momenta and translational energies. The source assumes small
relative drift and homogeneous fixed volume; it is not a moving field solver.
Initial argon velocities are +1000 m/s axial for neutrals, -1000 m/s axial for
ions, and 10000 m/s azimuthal for electrons, with heavy temperature 0.1 eV.
All charge states retain their physical neutral-minus-bound-electron masses.

Three complete 3.5 us protocols (all Ar+, 10% seed, heated 10% seed) run at
5/2.5 ns. The heated case has the existing 20 eV/nucleus pulse. A common
[20000,-10000,5000] m/s boost of all species is a frame-invariance comparison:
populations, temperatures and quasineutral current should be unchanged. It does
not certify accuracy of the Maxwellian or small-drift closure.

Independent conversion limits use declared constant reaction coefficients,
not fitted argon data. Ionization has x(t)=1/[1+(1/x0-1)exp(-N*k*t)],
k=1e-14 m³/s. Recombination has x(t)=x0/(1+N*alpha*x0*t), alpha=1e-18 m³/s;
its constant Maxwellian rate implies mean capture energy 3*Te/2. Analytic
per-event donor/product momenta and two-body centre-of-mass energies determine
all expected species moments. Full 50/100/200-step histories measure convergence
against those independent formulas. No assertions derived from the generator.

A moving excited-neutral single-line decay has exactly zero electron number,
electron momentum and electron energy. Population, transferred heavy momentum
and energy follow exp(-A*t); photons carry the independent binding release.
This exposes any fabricated seed or lost motion of initially empty states.

Outputs retain all accepted states, measured number/charge/mass/momentum/energy
budgets, currents, physical ordering ratios and the source archive. Rejections
remain failures, without floors or repairs. Private OPEN-ADAS input paths and
hashes identify data; provider tables and rate matrices are never bundled.
Missing channels/branching and photon absorption, spatial/current feedback,
neutral scattering and finite-drift reaction physics remain full-goal work.

The first run outputs/atomic_moments preserves actual underflow failures and the
first rejected endpoint arrays. Multiplying representable tiny populations by
kilogram masses caused a zero denominator. Its postprocessing also failed when
comparing unequal rejected prefixes. A separate process_exit.json records the
actual terminal exit; revised comparisons skip failed full trajectories and
retain their failures. Fresh runs never overwrite this evidence.
