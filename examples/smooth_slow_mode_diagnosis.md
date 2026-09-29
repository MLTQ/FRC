# Slow-mode magnetic cancellation diagnosis

Linearizes about the original uniform n and axial B0. Independently averages
true linear bulk momenta and analytic perturbation-only psi/g, then varies only
the velocity/current used in -v×B0+delta J×B0/(ne). Four controls use numerical
velocity/current, one exact quadrature field at a time, or both exact. The
uniform-density pressure gradient is curl-free and is explicitly omitted only
from these magnetic controls. Physical finite inertia still acts in every solve.

The perturbation-only weak current avoids bias subtraction noise in this linear
diagnostic; production states are never reinitialized or altered. Ideal/Hall-to-net
field ratios measure cancellation independently of the numerical derivative.
Complete rate arrays and errors are retained. Exact substitutions are diagnostic
oracles, not acceptable formation-solver patches or claimed full trajectories.
These controls distinguish magnetic-space representation and velocity/current
mapping errors from the corrected conservative thermal transport.
Reference metadata uses only declared dataclass parameters, excluding cached
eigenvectors/intermediate arrays populated during evaluation.
