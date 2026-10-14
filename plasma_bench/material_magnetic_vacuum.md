# Magnetic vacuum following a material boundary

`MaterialMagneticVacuum` extends the finite-cylinder exterior to actual free
side/end motion. Only magnetic flux exists outside the plasma: no mass,
conductivity, pressure or time-evolved vacuum current is added. The fixed outer
rectangle retains the existing continuous multipole trace and infinite energy.

Boundary flux coefficients use the parent's polynomial **psi** trace, regular
at the axis. Future plasma coupling must use that same trace. The existing
compact vector space interpolates r0²*a instead, and must not be spliced into
this exterior at the ends merely because nodal values agree. This module alone
does not complete moving plasma/current/coil coupling.

Q1 material boundary displacement extends linearly through each exterior patch
to zero at the outer box; corner patches use the product of both blends. This
is a computational coordinate rule, not vacuum inertia. Folded geometry rejects
without mesh repair. All side and end coordinates, including axial pole motion,
participate. Rotation does not change the axisymmetric poloidal vacuum geometry.

For meridional F, actual r and physical flux gradient g=F^-T grad0(psi),

    H = (2pi/mu0) integral det(F)/r [g_r.T*g_r+g_z.T*g_z] dR dZ + H_infinity.

Vacuum internal nodes/multipoles minimize this energy at the actual boundary
flux. `geometry` returns the condensed Hessian and extension. At specified
boundary velocity, its exact geometry derivative uses

    metric = (div2(v)-vr/r)*I - grad2(v) - grad2(v).T.

The minimized Hessian rate is extension.T*Hdot_full*extension; derivative of
the minimizing interior extension cancels by stationarity. `evaluate` supplies
physical Br/Bz, volume, energy, boundary work covector and coordinate energy
gradient. `field_response` reuses an already assembled geometry during a coupled
plasma stage, avoiding a second vacuum elimination. The derivative is the analytic
derivative of the SAME geometry-dependent
energy, pulled back through the mesh-motion map. Its negative is the exterior
force contribution at held boundary flux, not a residual-fitted force.

The finite-space vacuum may have shape sensitivity to the interior coordinate
rule. Independent analytic fields, exterior energy and Maxwell traction on a
moving boundary must measure its spatial convergence. Exact algebraic work
consistency does not establish physical accuracy. No temporal/spatial smoothing,
diffusion, field substitution or energy repair is applied. Fixed-cylinder
`samples` is rejected to prevent accidental reconstruction on obsolete geometry.
