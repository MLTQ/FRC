# Complete low-level argon radiative trajectories

`python -m examples.line_radiation NEW_OUTPUT` runs complete spontaneous-only
cascades using the 15-level, 31-line NIST record. Use the established headless
plot environment. Raw input, source closure, every accepted state, line photon
counts, plots and numerical measurements are preserved with existing tools.

`evolve` uses RK4 for level fractions and cumulative photon counts per line.
It never corrects population sum or energy. `sequential_reference` independently
evaluates the exact Bateman solution of the physical 2p5→1s4→ground cascade,
using the published 4e7/1.32e8 s^-1 probabilities. Complete 50 ns trajectories
at 1/.5/.25 ns compare all three populated levels, not just a final invariant.

Three complete 3.5 us cascades at the same steps begin with 90% ground-state
atoms and 1% in each of the ten 4p levels. This is a declared prepared
distribution to isolate radiation, not a prediction of electron excitation or
laser/plasma preparation. Four 4s levels evolve separately; finite-lived
resonance states radiate and the truncated model retains metastables. Physical
metastable loss channels remain to be supplied. Photon production is measured
separately from stored energy, and is not called escaped energy.

The same atomic data evaluate line-centre absorption for uniform, ground-state
neutral columns at 300 K and 0.1–100 mTorr, using a 7.5 cm centre-to-side path.
The JSON also gives 0.5 m centre-to-end columns at the prior source-study
nucleus density, 0.1–100% ground-neutral fractions, and 300 K/1 eV atom
temperatures. Pressure equivalent is a neutral density label, not the total
plasma pressure. Profiles include thermal plus natural widths only.

These optical depths constrain the optically thin assumption; they are not
escape probabilities, radiation transport, collisional rate measurements,
gas-fill recommendations or FRC formation results. Opacity is not used to
multiply A by an arbitrary escape factor. The moving pulse closure is unchanged.
