# Full original spatial RHS contribution records

## Purpose
Decompose the original tangent probe's twelve accepted signed RHS evaluations
without changing a physical state, source equation, wall space or solver bound.
The design and equations are in docs/SPATIAL_PARTS_DESIGN.md.

## Components
- `species_parts`: reproduce each species' exact affine pressure split and
  compare their sum bitwise against the unchanged production pair.
- `decompose`: archive ten fluid/forcing groups, solve the eight induction
  groups, retain rejections, full species-restored Ohm checks and wall reactions.
- `run`: hash-check original source/raw replay, reproduce full rates bitwise,
  perform all signed cases and save complete term/forcing/input records.

## Contracts
The original fixed1e-10 induction tolerance is required for every component.
An individual rejection is recorded, with no manufactured missing rate.
Fully returned sums retain norm differences and original Ohm/Ampere/Faraday
backward errors. Full matching does not prove spatial continuum accuracy.
Species forces are derived directly from dh,ds, not differences of large fields.
Actual passed inputs are checked unchanged. Both native wall reactions remain
explicit. Component coordinate projections omit inactive wall-normal slots
because the existing declared coordinate space does; physical audits do not.
Full physical comparisons and a fresh complete workflow replay are required.
Used source becomes frozen. No old source, test, time integration or state changes.
