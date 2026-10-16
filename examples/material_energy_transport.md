# Complete numerical heat-transport limits

Runs the same conservative thermal-energy diffusion function used in the moving
pulse against independent analytic heat equations, then preserves all histories,
source archive, measurements and a scientific PNG. This is a concrete numerical
transport comparison, not a unit-test suite, reporting framework or claim that
this numerical diffusion supplies the missing physical heat conductivity.

Periodic unit-length domain with initial cell averages of h=1+0.8 cos(2*pi*x).
For constant D, the independent continuum solution multiplies the sinusoid by
exp(-4*pi^2*D*t); exact cell averages include sinc(dx). Fixed D=0.1 gives a
second-order spatial diffusion comparison at16/32/64/128 cells. Time evolution
uses the actual graph exponential through51 saved states to t=0.5.

A second set uses fixed species drift |u|=1 and D_num=|u|*dx/2, the method's
declared upwind-scale numerical coefficient. Its continuum comparison remains
explicitly diffusive at each grid. Distance from the initial profile measures
how that artificial transport vanishes as spacing tends to zero. These cases
evolve the diffusion component alone; they do not silently omit pressure/current
feedback from a purported coupled plasma experiment.

An unequal-particle two-reservoir case has m=(1,3), h0=(4,0.2), g=0.7.
Its independent solution is h(t)=h_mean+(h0-h_mean)*exp(-g*(1/m0+1/m1)*t).
All cases retain total energy and physical entropy observations without clipping,
normalization or fitted rates. Whole moving-pulse time/space comparisons remain
necessary to assess central transport, split coupling, shocks and the edge.

Completed run93705 (0.891s), source84047d488b6f16f9c86bd37aeea54653a3d3dcedddb5400319218568dec8f3a8:
fixed-D final L2 errors at16/32/64/128 cells are1.99547e-3,4.98472e-4,
1.24592e-4,3.11464e-5 (second-order decrease). Maximum measured relative energy
change is2.73115e-14; every sampled total-entropy increment is positive.
For fixed drift the change from the initial profile decreases0.25635,0.14948,
0.08075,0.04197 as the numerical diffusion vanishes. Unequal-reservoir maximum
absolute solution error is5.32907e-15 and energy change1.55431e-15. Full histories
and the inspected figure are in outputs/material_energy_transport in the main
checkout. These results concern the added mixing term, not the whole plasma.
