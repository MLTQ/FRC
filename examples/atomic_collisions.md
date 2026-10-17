# Complete reacting/colliding species histories

First advances the physical electron/argon collision-only problem already solved
in maxwellian_pair_relaxation. Its archived independent velocity-space quadrature
solutions at 100/200/400 steps are reused without recomputing a reference or
building a new verification framework. Atomic conversions are exactly disabled
in this declared collision limit. The original ion mass (39.948 u), density,
large drift, temperature and common boost are retained; neutral mass is set to
ion mass plus electron mass solely to reproduce that existing collision case.
This does not change masses in the physical reacting argon case.

Then the complete 3.5 µs heated 10%-seed argon protocol runs at 5/2.5 ns with
reactions and all occupied charged-level/electron collision pairs together.
The physical argon neutral mass is 39.948 u, charged masses subtract q me,
and lnΛ=10 is explicitly prescribed. Actual level velocities/temperatures enter
every pair, rather than assuming one common ion fluid. The default omits neutral
elastic scattering and charge exchange; the explicit option below adds only
electron/ground-neutral elastic scattering.

All accepted extensive states, currents, species heat, independent budgets and
the original source identity are retained. The collision-only limits compare
full time histories to the existing independent integrals. The argon pair reports
time-halving differences and a physical comparison against the completed
reaction-only history. Failed prefixes are labeled incomplete, never clipped or
compared as though complete. No unit tests or exhaustive replay layer is added.
Derived states and input paths/hashes are saved; provider tables are not copied.

This evolves local current and heating self-consistently with changing charged
populations. Spatial transport, coil/field feedback, remaining neutral channels,
kinetic validity, radiation transport and real boundaries are still required.

The first unsplit attempt is preserved in outputs/atomic_collisions, including
the actual first rejected state with negative rare-charge thermal energy. The
current full step uses symmetric collision-half / reaction-full / collision-half
composition, retaining every pair and the same physical reaction inputs. Its
collision-only and full reacting histories are rerun in a fresh output directory.

`--neutral-elastic-file PERSONAL_LXCAT_TEXT --charged-reference COMPLETE_RUN`
adds the independently evaluated electron/ground-Ar elastic channel to the
original full3.5us heated protocol, using `ElasticAtomicMoments`. Ground targets
alone use the provided cross section; excited-neutral populations without a
matching channel are measured, not silently assigned the ground rate. The report
records their maximum fraction, input identity/domain and maximum observed
relative-velocity probability outside the table. That probability does not bound
unknown tail force/heat. Every retained charged pair and atomic channel remains.

This option reuses the already completed charged-only limits without rerunning
them; it executes the full700/1400-step heated histories, compares their common
times, and compares the fine result with the complete charged-only and reaction-
only histories. The charged-only reference must be terminal and successful.
New numerical source evidence is separately archived; no provider table is copied.
The comparison also includes array-assembled reactions whose full-history
equivalence was established separately. Source splitting, distribution closure,
unknown neutral channels and incomplete radiation physics remain limitations.

`plot_histories` renders saved histories independently of evolution. Current
magnitude and both nonzero vector components are shown: faster azimuthal damping
can coexist with a persisting axial current driven by the different neutral/ion
flows. A shared outside legend preserves the data panels. Replotting completed
states does not rerun them; the simulation and figure source identities remain
separate when the renderer changes after a completed experiment.
