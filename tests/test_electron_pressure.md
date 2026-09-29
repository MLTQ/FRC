# Electron charge, pressure and thermal tests

Independent physical current, direct analytic Biermann/entropy derivatives,
rigid electron heat translation and adaptive DOP853 evolution verify the new
split component. Checks include axis/current boundary traces, local charge
cancellation, actual curl of gradients, work with nonzero wall Poynting input,
and agreement with the prior Hall adjoint on its restricted electric space.
The gradient-free conditioned pressure curl is checked against the full electric
curl and independently differentiated continuum fields.
The omitted-wall-trace control checks the independently known surface-current
offset and wrong sign, rather than an arbitrary mesh-dependent error threshold.

Actual drift propagation refines on three meshes; the magnetic field remains
fixed and electron pressure moves with ve while bulk flow is absent. Time order
is checked against a separate adaptive ODE, not another midpoint step. Energy,
linear entropy proxy and physical logarithmic entropy are checked separately,
including the finite-step thermodynamic discrete-gradient identity. Logarithmic
means are independently compared to high-precision Decimal calculations.
Work cancellation uses its actual absolute summation scale; a nearly vanishing
net exchange is not an appropriate floating-point backward-error denominator.
Adiabatic and isothermal barotropes and proxy density-scale changes are controls.
Invalid boundary closures and failed solves preserve all input physical arrays.
