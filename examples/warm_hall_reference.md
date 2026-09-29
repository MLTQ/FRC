# Independent warm Hall linear continuum oracle

The background density and both species pressures are uniform. At first order,
div J1=0 implies div ve1=div vi1. The adiabatic species equations therefore give
pe1=fe p1, pi1=(1-fe)p1 with cs²=gamma(pi0+pe0)/rho0. The existing independently
derived finite-sound-speed Hall generator/eigenvalues consequently apply with
**total** pressure. This is a derivation under equal gamma/uniform background,
not an assertion that arbitrary warm Hall induction is unchanged.

`WarmHallWave` includes the physical electron velocity vi1-J1/(n0 e) and adds
the exact Bessel gradient -grad pe1/(n0 e) to Ohm's law. Its curl is zero only
at linear order about this uniform background. Species enthalpy drift has
zero divergence at this order but remains essential in the nonlinear solver.
The summed quadratic species-pressure energy equals p1²/(2 gamma p0).

J1(qR)=0 sets Br/Btheta/vr/vtheta=0 at the radial wall; psi is fixed there.
The warm pressure Ez(R) is generally nonzero, and is allowed by the explicit
zero-Btheta magnetic trace. These fields are not full conducting-wall electric
data. Periodic z, fully filled cylinder, no collisions/resistivity/inertia.
Finite-amplitude fields have quadratic nonlinear residuals; amplitude controls
are required for actual evolution error. fe=0 is the exact cold reference limit.
