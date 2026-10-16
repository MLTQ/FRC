# Full independent affine limit with explicit source subcycling

Repeats the unchanged16-cell warm neutral expansion/decay physical problem
through its full6.083704µs endpoint, now with20 explicitly declared reaction
intervals. Each interval retains fine actual acoustic transport steps. Mass,
temperature, radiation rate, precise geometry and initial fine density are the
same as examples/precise_material.md. Physical-width coarsening is zero.

The independent affine energy integral fixes front/velocity at scale1.5;
population and photon energy follow spontaneous exponential decay. This neutral
limit has commuting common motion/source physics, so agreement cannot certify
the source interval in the charged case. It does directly exercise all source
halves, internal transport steps, variable support return and reservoir budgets
through a complete trajectory. No unit tests or custom verification framework.

Outputs retain full completed coupled states, all internal dt/budget records,
source archive and measured analytic errors. A failed interval remains failed
with the original duration requested. Use `python -m examples.subcycled_affine
FRESH_OUTPUT`. Charged source-interval refinement remains separately required.
