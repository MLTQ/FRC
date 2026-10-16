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

## Complete unboosted histories and preserved frame failure

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

Revised run 7072 exited 1 in 832.137 s in `outputs/atomic_moments_reordered`, source
5453c1ea52380e61dd5330e4cc5947df9eb65f036f5fa85eaa9eb1fb03170174.
All six conversion limits, all six unboosted argon histories and the exact-zero-
electron moving neutral decay completed. Only the commonly boosted heated case
failed: 79 accepted steps, rejection at attempted step 80. Its actual rejected
state and explanation are in boost_rejected_state.npz and
boost_underflow_diagnosis.json. Heavy population 69 and its momentum round to
zero while its lab energy is 5e-324 eV/nucleus. Multiplying by a large common
bulk energy exposes the finite-range limit; no value is rounded away to pass.

| Protocol, 1400 full steps | Final electrons/nucleus | Final Te, eV | Maximum energy drift, eV/nucleus |
|---|---:|---:|---:|
| Initially all Ar+ | 1.049276875 | 2.317779959 | 3.93641e-12 |
| 10% seed | 0.135647194 | 2.386806768 | 4.38760e-13 |
| Heated 10% seed | 0.683902917 | 2.536280756 | 2.03260e-12 |

Across the 5/2.5 ns pairs, maximum population-L1 differences are
7.46246e-6 / 2.15925e-6 / 1.72802e-5 for these three protocols. Corresponding
maximum Te differences are 0.000138737 / 0.000182440 / 0.001049801 eV, and
current-vector differences are 0.000109936 / 8.01661e-7 / 1.85767e-6 A/m².
The small-drift ordering reaches at most 0.0106079 across these histories.
The neutral decay has zero electrons and zero current throughout, with maximum
independent population/momentum/energy errors 1.88738e-15 / 4.44089e-15 /
1.99840e-15 in the stated coordinates. The scientific PNG was inspected.

In the complete constant-ionization reference, current changes by at most
3.82e-11 A/m² while electron number rises from0.2 to0.648783. This has a
physical interpretation: each event adds equal and opposite electron/ion charge
flux at the donor velocity. Carrier growth reduces relative kinetic energy
through the evolved particle moments, rather than by an added energy correction.
The separate recombination reference removes current as charged particles merge;
its final electron temperature remains2 eV under the constant-rate capture law.
All accepted histories and explicit failures are preserved; incomplete boosted
prefixes are not presented as a complete frame comparison. The separate
`FixedFrameAtomicMoments` now offers fixed translating coordinates, with exact
nonrelativistic moment transformations and laboratory current/energy output.
Read plasma_bench/atomic_frame.md. The original source and its raw failure remain
unchanged. Run 95433 exited 0 with both full boosted heated trajectories in
distinct fixed frames, preserved in outputs/atomic_moments_fixed_frame under
source 8ecabaefe7b52d082521657f438dfa851a6552732d71aed5ec8316abaa3d7734.
The first frame equals the common boost [20000,-10000,5000] m/s; its populations,
Te and transformed current agree exactly with the unboosted trajectory because
the integration-frame initial states coincide. This alone would be a weak check.
The second frame adds [137,-211,93] m/s and uses different extensive states.
Over all 1400 steps its maximum differences from the transformed reference are:
population L1 1.25376e-14, Te 5.32907e-14 eV, current 2.25265e-9 A/m²,
species momentum 1.02141e-13 in the stored scale, and lab energy 8.59757e-13
eV/nucleus. These are coordinate-consistency results, not new atomic accuracy
or timestep-convergence evidence. Both complete histories and the inspected
comparison PNG/PDF are preserved. Finite float range and negative-thermal-state
rejection remain explicit; choosing a frame does not cure arbitrary relative
velocities or remove physical model limitations.
Do not install this source into the moving pulse based only on these homogeneous
balances: spatial transport, field feedback and large-drift validity remain open.

OPEN-ADAS raw tables remain personal external inputs with recorded paths/hashes;
there is no provider-data bundling, LXCat terms acceptance or publication.
Bead frc-transport.5.5 tracks both these moment equations and the remaining
physical/spatial/current coupling. The overarching FRC goal remains open.
