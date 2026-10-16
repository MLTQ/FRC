# Complete material-fluid wall-shock trajectories

Independent physical baseline for the optional numerical bulk viscosity. Uniform
Ar-mass gamma5/3 gas flows at10km/s into a fixed axial wall; p0/(rho0*u0²)=0.01.
Radial motion is constrained by rigid sides, the entering left boundary moves at
u0, and its constraint reaction work is observed. This is explicitly a bounded
one-dimensional fluid fixture, never a plasma formation experiment or vacuum.

The same consistent material mass, cylindrical volume, bulk-viscosity force and
dual-species heating are used. Restricted mass solves enforce boundary velocities
without momentum/energy repair. Numerical heat goes to ions. No electromagnetic
or exchange terms are included because this is the neutral Euler shock baseline.
Four radial cells satisfy the shared mesh contract; radial motion is constrained
at every node. Each of50/100/200 axial cells evolves all400/800/1600 RK4 steps to20us. All final physical
arrays, independent work and volume-weighted L1 differences are saved.

The exact shock is derived from mass, momentum and energy jumps with downstream
velocity zero: a=(gamma-1)u0/2, c0²=gamma*p0/rho0,
s=[-(u0-a)+sqrt((u0+a)²+4c0²)]/2. Front=z_wall-s*t, rho2/rho0=(u0+s)/s,
p2=p0+rho0*u0*(u0+s). It does not use numerical solver output or fitted parameters.
The startup boundary layer and wall-heating error remain visible. This is a
complete simulation comparison, not a unit-test or general verification framework.

`--cells 100 --time-refinement 4` repeats the complete fixed-grid physical shock
with four times as many timesteps. It separates the entropy-heating time error
from the spatial shock width; it never changes or fits the analytic comparison.
