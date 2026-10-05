# Axis and wall local-budget norms

Uses complete hash-verified physical cell-budget arrays and independently
reconstructs the current-gradient normalization at each quadrature. Reports
relative/absolute L2 and maximum cell-density residual in the first and last
fluid radial elements, for particles, current, both entropies and total energy.
No local residual is dropped because of small annular volume or global
cancellation. The first/last masks follow declared aligned mesh topology.

Preserves source/input identities and can replay every metric bitwise. This
postprocessing leaves the original frozen flux calculation untouched. It is a
regional cell-integral audit, not a claim of pointwise conservation between
quadrature points or a physical repair of current-normal jumps.
