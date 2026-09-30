# Independent time-dependent cylindrical Euler manufacture

Density, both pressures and all three velocities vary in radius, axial position
and time. The radial velocity is odd and zero at the outer wall. Scalars are
even and angular momentum vanishes quadratically at the axis. The pressure
gradient at the wall balances centrifugal force at every z/time, so forcing
does not conceal an incompatible radial wall acceleration.

`source` evaluates the **continuous** conservative cylindrical Euler residual
using complex spatial/time derivatives, independently of production derivatives
or fluxes. Its radial geometric source is (p/M+n vtheta²)/r; angular momentum
has no geometric torque. `primitive_source` independently derives the same
residual from original material velocity and species adiabatic-pressure laws.
Pressure is differentiated directly; the analytic cancellation of p/r is
performed before floating-point evaluation, avoiding a weak radial source
being lost between two large geometric pressure terms near the axis.
Both evaluate only positive volume radii; exact axis regularity follows from
the analytic parity, not a small-radius floor.

The imposed source exists only for verification. This exact time-dependent
field is not a free Euler solution, an FRC formation experiment or a physical
external heating/force model. The unforced companion experiment uses its
initial state but has **no analytic time-dependent target**.
