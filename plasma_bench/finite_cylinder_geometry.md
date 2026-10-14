# Full finite-cylinder compatible electromagnetic space

Poloidal flux has independent nodes at both physical axial ends; only the regular
axis is fixed to zero. Er/Etheta use axial nodal polynomials; Ez and Btheta use
axial edge histopolants. The full Faraday incidence therefore has finite end
traces instead of periodic wrapping. All current/electric physical masses use
the actual cylindrical quadrature and full theta boundary values.

Vacuum energy is the FiniteCylinderVacuum Hessian mapped from exactly the same
polynomial side/end flux trace. Its work transpose enters the same full magnetic
Hessian as interior energy. No exterior matter/current/conductivity appears.
Toroidal vacuum B is zero: the interior toroidal space has explicit radial-wall
and both-end trace null spaces. Its declared Galerkin Faraday projection retains
curl(gradient)=0 and does not repair accepted states.

The magnetic field uses the same regular flux convention as before. Pack/unpack
no longer duplicate an axial endpoint. Axis regularity, finite arrays and all
three toroidal boundary traces are checked. Old periodic geometry remains a
separate control. Finite current/fluid maps and the coupled integrator must use
these changed dimensions; old periodic coordinate maps cannot be reused.

Seven full24us coil pulses now exercise this coupled path. Independent limits,
spatial/time comparisons, transport flags and remaining exterior-resolution
limits are recorded in docs/FINITE_CYLINDER_RESULTS.md.
