This separate variant tests the radial thermal/kinetic work adjoint. It retains
the complete original protocol and records, with its own source closure. The
failed baseline and its exact audit source remain unchanged.

# Self-consistent fluid trajectory protocol

Fourteen unforced cases: contact/full transit, pre-shock finite-amplitude simple
wave and solid rotation at nr8/16/32, plus nr16 half-step, reverse-flow and
higher-quadrature controls. z has twice the radial cell count. The fixed RK4
step is chosen from both mesh dimensions and 1.5 reference sound speed; this is
a declared fixture choice, not a universal CFL/positivity theorem.

Keeps every accepted global history, five complete snapshots, initial rates,
last pre-step state and local integrated spatial energy defect. Actual physical
energy and conservative energy-account drift are distinct. Failed steps retain
the actual last accepted state/history and traceback without repair or retry.
Initial data are cell averages and nodal entropies of independent exact fields.
Initial projection and subsequent error are both recorded. No analytic forcing
is applied during evolution. Norms use an independent physical sampler.

`--cases` selects indices without overwriting an existing output directory.
Two worker processes by default; source closure is frozen and rechecked before
completion. The protocol is a bounded fluid verification, not FRC formation.
