# Moving resistive verification measurements

Records 20 independently initialized moving cases: FRC/source mesh refinement,
a wider physical diffusion layer, outer time refinement, matching-sphere,
diffusion-tolerance, mode and quadrature controls. A different scalar diffusivity
changes physical equations; it is not tuned numerical stabilization.

Accepted histories retain cumulative physical heat, source work, full response
error and adaptive diffusion counts. Initial/mid/final lossless arrays include
evolving entropy, flux labels and actual current. Partial/failure records stay
explicitly incomplete and are never assembled as successful studies. Reuse
requires matching full source/configuration and completed data. Disjoint named
cases may run in separate processes; no two jobs should own the same output.
