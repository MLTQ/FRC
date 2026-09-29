# Independent spherical mode with electron inertia

The existing analytic sphere mode keeps its continuous dipole/vacuum matching.
With de²=me/(mu0 ne e²), its decay rate becomes D k²/(1+de² k²).
This follows by substituting Jtheta=k² r a/mu0 into the continuum Ohm law,
independently of the finite-element matrices. Magnetic and current kinetic
energies have ratio 1 : de² k²; integrated Joule heat loses their sum.

`InertialSphere` reuses only the independent analytic spherical radial function.
Physical diffusivity and effective decay diffusivity are distinct. Local heat
and total heat use physical eta, not the smaller effective decay coefficient.
The axis, all-space energy and interface contracts are inherited unchanged.
Zero electron length explicitly selects the massless reference.
The underlying analytic radius, diffusivity and mode contract validates inputs
before computing the effective rate.

The example is frozen geometry, linear toroidal induction, not mechanical force
balance. Small field amplitudes make nonlinear motion/Hall terms higher order
on the comparison interval; this does not validate finite-amplitude formation.
