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
every pair, rather than assuming one common ion fluid. Neutral elastic scattering
and charge exchange remain absent, not established negligible.

All accepted extensive states, currents, species heat, independent budgets and
the original source identity are retained. The collision-only limits compare
full time histories to the existing independent integrals. The argon pair reports
time-halving differences and a physical comparison against the completed
reaction-only history. Failed prefixes are labeled incomplete, never clipped or
compared as though complete. No unit tests or exhaustive replay layer is added.
Derived states and input paths/hashes are saved; provider tables are not copied.

This evolves local current and heating self-consistently with changing charged
populations. Spatial transport, coil/field feedback, neutral collision data,
kinetic validity, radiation transport and real boundaries are still required.

The first unsplit attempt is preserved in outputs/atomic_collisions, including
the actual first rejected state with negative rare-charge thermal energy. The
current full step uses symmetric collision-half / reaction-full / collision-half
composition, retaining every pair and the same physical reaction inputs. Its
collision-only and full reacting histories are rerun in a fresh output directory.
