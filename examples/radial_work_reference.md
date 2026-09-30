# Two-dimensional original Euler rate reference

Independent analytic density, both pressures and all three velocity components
vary in r/z. Complex-step spatial derivatives inserted into original cylindrical
Euler/adiabatic pressure equations supply initial physical rates. Normal wall
velocity is zero and the pressure derivative there balances centrifugal force.
This is an instantaneous rate reference, not an analytic evolved solution or
a forced trajectory. No production reconstruction or flux supplies its target.
`derivatives` supplies the same spatial derivatives plus those PDE time rates
for the common physical norm calculator. Temperature and energy rates follow
the independently differentiated physical EOS and kinetic energy.
