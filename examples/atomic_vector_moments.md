# Complete reaction histories with unchanged resolution

Runs the array-assembled moment source through the existing independent full
constant-rate ionization (2 µs) and recombination (20 ms) limits at 50/100/200
steps. Expected populations follow the logistic/reciprocal laws; donor and
captured-pair momenta/energies follow the independent per-event formulas in
atomic_moments.md. The same existing `evolve` routine records all actual states,
physical budgets and full-duration failures. No unit suite is used.

Four real-rate ground-ladder argon histories compare the old and new assembly
at 700/1400 steps through the full 3.5 µs heated protocol. Initial 10% seed,
10 eV electrons, 0.1 eV heavy species, neutral/ion ±1000 m/s and electron
10000 m/s transverse velocities are unchanged. The prescribed pulse deposits
20 eV per initial nucleus. This is the public ground-ladder source used in the
moving argon attempt; it does not add unresolved excitation or private tables.

Each complete old/new pair reports elapsed wall time and direct all-state
population, momentum, energy and current differences. These are numerical
equivalence/performance observations, not experimental validation or evidence
that the underlying closure is adequate for formation. Timings include the
existing physical observations and depend on concurrent machine activity.
The independent conversion comparisons provide separate mathematical accuracy.
Raw histories, input hashes, source archive and a scientific PNG/PDF are saved.

Usage: `python -m examples.atomic_vector_moments FRESH_OUTPUT`.
