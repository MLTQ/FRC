# Collision sensitivity of an actual completed pulse

Evaluate the finite-drift paired Maxwellian collision moments at each recorded
state of a completed configured experiment. The scalar Coulomb logarithm is
explicitly 10; its screening/cutoff prescription has not been established for
this pulse. Compare original eta*J², small-drift Maxwellian drag using the same
logarithm, finite-drift drag and separate electron/ion heat. No field, temperature
or current is changed, and no candidate trajectory is evolved.

The integral along this **original** trajectory is a sensitivity diagnostic;
it is not predicted heating of a future pulse using altered transport. Its sparse
summary-state quadrature is also distinct from the original solver's every-step
Joule integral. Each sample's heat sum and drag are retained without correction.
Maximum local relative-drift relaxation rate times original dt screens one
source timescale; it is not a full-system stability criterion or timestep choice.

The example freezes its own source/companions and hashes input states/summary,
then emits one PNG/PDF plus actual values. The model-validity limitations in
`maxwellian_pair_collisions.md` remain. This is physical model exploration, not
a pass/fail suite, convergence certificate or new verification framework.
