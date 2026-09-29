# Velocity reconstruction at compatible physical quadrature

Fits actual conservative cell averages with degree-3/5 annular polynomials
r^s P(x), x=r²/2, where s=0 for rho/mz, s=1 for mr and s=2 for angular
momentum. Two mirrored stencils are averaged, with physical interior clipping
at radial boundaries and periodic axial averaging. Density is reconstructed
before dividing momenta. No analytic field/current enters that reconstruction.

`velocity` returns physical (vr,vtheta,vz) at the supplied geometry's common
Gauss points. No continuity flux or stored density/thermal energy is changed.
Negative reconstructed density rejects; there are no floors/limiters. This is
an instantaneous smooth field map, not a full moving velocity/vorticity/energy
contract. Boundary nodal traces, nonlinear force/work, vorticity and coupling
still require consistent derivation before adopting it in the moving solver.
