# Atomic reactions with species momentum and energy

The homogeneous joint atomic network now has an optional moment model in
`plasma_bench/atomic_moments.py`. Each resolved heavy level and the electron
fluid carries its own momentum and total translational energy. Changing charge
therefore changes particle mass, current and kinetic energy explicitly. This
component is not yet connected to moving geometry or electromagnetic evolution.

The physical model and state equations are in the module companion. Ionization
uses donor-velocity moments under the small-drift Maxwellian approximation of
[Meier and Shumlak](https://www.aa.washington.edu/sites/aa/files/research/cpdlab/docs/Meier_PoP2012.pdf).
Radiative capture uses independently derived two-body centre-of-mass kinematics:
relative kinetic energy and binding release go to the produced photon reservoir;
the product retains the pair centre-of-mass motion and its random component.
Photon recoil/momentum corrections, finite-drift selection, neutral scattering,
charge exchange and the existing atomic-data/opacity limitations remain explicit.
The large drift of the strong moving pulse can violate this reaction closure.

Masses are neutral atomic mass minus q electron masses, with nonrelativistic
binding mass defects neglected. Separate moment coordinates for every heavy
level retain motion when an initially empty state is populated. Fixed-state
reaction generators conserve nucleus count, charge, physical total mass, vector
momentum and lab translational + atomic + produced-photon energy. Conservation
is observed, never restored by overwriting density, temperature or a reservoir.
The nonlinear exponential-midpoint integrator remains second order.

## Independent complete conversion limits

The first and reordered implementations both completed six mathematical
conversion trajectories. These deliberately use constant coefficients, not
argon fits: ionization k=1e-14 m³/s, recombination alpha=1e-18 m³/s, N=1e20 m^-3.
The independent number solutions are logistic ionization and reciprocal
recombination. Per-event donor/product kinematics independently determine all
expected species momenta and energies. No unit-test suite was used.

| Full steps | DI population error | DI lab-energy error, eV | RR population error | RR lab-energy error, eV |
|---:|---:|---:|---:|---:|
| 50 | 3.76923e-5 | 5.12612e-4 | 2.58900e-5 | 7.76726e-5 |
| 100 | 9.50287e-6 | 1.29238e-4 | 6.39717e-6 | 1.91921e-5 |
| 200 | 2.38586e-6 | 3.24475e-5 | 1.58973e-6 | 4.76934e-6 |

Finest maximum momentum-coordinate errors are 4.87710e-7 (DI) and 1.62395e-7
(RR), with momentum scaled by sqrt(m_neutral*eV). Errors decrease approximately
fourfold under halving. Across these complete limits, total mass error relative
to neutral mass is <=3.00e-15 and total momentum-coordinate drift <=2.77e-16.
These analytic limits verify the conditional source equations, not the full
argon kinetics, large-drift reaction physics or laboratory applicability.

## Preserved representability failure and current run

Initial run 94014 exited 1, source
7640a79e64bc32c0745d2de87543817b4df7459a3bf242c52d977b4b66d6830e.
`outputs/atomic_moments` preserves six complete analytic histories, two complete
seeded argon histories and a moving-neutral decay, plus failed all-ion/heated
prefixes. Postprocessing subsequently failed on unequal failed-prefix lengths;
`process_exit.json` records the actual process termination separately from the
unfinished original summary. Do not treat its nonterminal summary as a live job.

The first rejected all-ion endpoint is saved in first_rejected_endpoint.npz;
the preceding accepted midpoint is separate. A nonzero population 4.8e-322
multiplied by its kilogram mass underflowed to zero before division, producing
an invalid velocity/thermal energy. State arrays and diagnosis are retained;
no provider generator/rate arrays were saved. Reordered arithmetic divides
dimensionless momentum by population before applying the physical mass scale,
and computes kinetic energy per particle before multiplying by population.
No population, momentum or energy is clipped or replaced. Representability of
the complete remaining trajectories still needs their actual terminal evidence.

Revised run 7072 is LIVE in `outputs/atomic_moments_reordered`, source
5453c1ea52380e61dd5330e4cc5947df9eb65f036f5fa85eaa9eb1fb03170174.
Its six conversion limits completed; full argon/boost/decay runs are pending.
The previously failing all-ion 700-step argon trajectory has now completed:
final electron number 1.04927407/nucleus, Te 2.31778445 eV, maximum energy
budget drift 4.65405e-12 eV/nucleus, and maximum relative-drift/electron-thermal
speed ratio 0.0106079. This is one complete trajectory; the finer and other
protocols still require their actual results.

In the complete constant-ionization reference, current changes by at most
3.82e-11 A/m² while electron number rises from0.2 to0.648783. This has a
physical interpretation: each event adds equal and opposite electron/ion charge
flux at the donor velocity. Carrier growth reduces relative kinetic energy
through the evolved particle moments, rather than by an added energy correction.
The separate recombination reference removes current as charged particles merge;
its final electron temperature remains2 eV under the constant-rate capture law.
The example compares 5/2.5 ns argon trajectories, common-frame boosts, and a
moving excited neutral with exactly zero free electrons. All accepted histories
and explicit failures are preserved; failed full trajectories no longer trigger
misleading comparisons of unequal prefixes. Keep the imported atomic source
unchanged until this exact process terminates. Do not install this source into
the moving pulse based only on count/energy balance or the six simpler limits.

OPEN-ADAS raw tables remain personal external inputs with recorded paths/hashes;
there is no provider-data bundling, LXCat terms acceptance or publication.
Bead frc-transport.5.5 tracks both these moment equations and the remaining
physical/spatial/current coupling. The overarching FRC goal remains open.
