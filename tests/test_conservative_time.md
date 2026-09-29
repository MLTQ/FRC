# Coupled SSP temporal verification

Independent nonautonomous scalar equations with analytic exponential solutions
replace only the physical RHS. They test both energy and total magnetic flux,
while actual changing coil fields and magnetic-vacuum completion remain active.
Timestep halving must converge at the declared second or third order. These
manufactured source terms are integrator tests, not physical MHD experiments.

The same known power checks every cumulative budget's positive RK weights.
Accepted Euler-stage counts must match the selected method. An injected third
stage failure leaves state, time, work and limiter counters unchanged. Invalid
time orders are rejected explicitly. Physical shock/source studies separately
test this integration method with the actual MHD equations.
