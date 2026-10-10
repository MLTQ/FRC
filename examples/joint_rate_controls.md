# Full-workflow spatial controls

Loads terminal full signed workflows, validates raw/source identities and
enforces each case's original declared induction tolerance. Reports every
volume/axis/wall physical norm difference between quadrature16/20, retaining
absolute scales for tiny errors. Interrupted progress cannot be accepted.
An explicit failed solve remains failed and cannot produce a paired result.

Optional same-grid/reference pairs compare centered full fluid, magnetic,
current and electric coefficients across amplitudes or residual bounds.
Their Euclidean coefficient differences are explicitly separate from physical
error norms. Saves complete difference arrays and every input hash. Physical
error-norm changes alone are not a norm of the difference field.

Fresh complete replay of this postprocessor checks arrays and non-runtime
metrics. It neither reruns unit tests nor changes a solver/state/bound. Actual
paired runs must already be complete, and their own full replays remain
separate requirements. No new threshold, formation or stability claim.
