# Warm boundary and local energy verification

Full scalar gradients have zero curl and the free-wall-Ez current work matches
the independent dual current at every DOF. Constant-weight RT0 remapping retains
local divergence-free charge, including zero normal wall flux. Varying wall
density with zero Btheta permits nonzero pressure Ez but no boundary work.

An independent J1-zero cylindrical diffusion mode verifies the new magnetic
trace in an actual implicit step, along with positive heat and total flux work.
The evolved-field asymptotic check uses nr16/32/64. A separate measurement
retains the nonconverging instantaneous first-annulus diffusion derivative;
this evolution test does not assert pointwise consistency of that derivative.
Analytically differentiated nonlinear electron fields give the expected local
thermal-plus-magnetic rate. Conservative pressure Poynting/enthalpy divergence
and its local ion-energy residual both must refine toward that reference.
