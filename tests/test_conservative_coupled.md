# Conservative plasma/vacuum verification

Checks stationary uniform states under both magnetic boundaries, independent
harmonic-polynomial exterior traces, positive interface diffusion with O(h)
scaling and exact zero vacuum resistance. Independent advection/diffusion
solutions check smooth convergence and inactive limiting. The original strong
shock requires decreasing B error, greatly reduced frozen-flux defect, no false
interior reversal, exact fluid energy-flux accounting and decreasing independent
response-energy error. A deliberately disabled-interface/limiter case verifies
atomic failure. The prescribed-coil test measures response error relative to
actual source work, not just large initial gas energy. Larger recorded cases
exercise the complete-update limiter and retain disabled failures.
An independent manufactured instantaneous Poynting identity with uniform B
checks that interface electric fields also transport energy on adjacent interior
faces; zero initial current must not create a spurious thermal-energy derivative.
Subtracting the zero-electric baseline isolates that manufactured nonideal term
from roundoff in differentiating the quadratic uniform-field flux.
The original nr=96 shock activates local complete-update limiting; its measured
affected region must remain small, exclude the central independent-reference
window, and retain frozen B/rho, shared-fluid energy, mass, CT and vacuum
constraints. Failure checks include unchanged cumulative limiter maps. Larger
response-budget and fine-grid accuracy studies remain recorded measurements.
