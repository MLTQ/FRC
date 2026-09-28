# Declared numerical shock diffusion

`shock_resistance` computes eta_num=mu0 C h min(h max(-div(v),0), |v|+c_bound),
where h=min(dr,dz), c_bound²=gamma p/rho+B²/(mu0 rho), and C is explicit.
It returns the numerical resistance in the existing nodal work equations and
the maximum eta_num. Positive neighboring-cell averages map it to plasma nodes;
vacuum and exterior modal rows remain exactly zero. Axis divergence uses regular
cell centroids; constant linear radial flow gives the correct 2 v_r/r.

The coefficient is **numerical**, not physical conductivity/transport. It is
O(h²) in smooth compression, O(h) across a compressive shock and zero in pure
expansion or resting fields. Positive resistive work is deposited as separately
reported numerical magnetic heating. This dissipative regularization addresses
the centered magnetic subsystem's shock oscillations; it is not a cell floor,
energy repair, or a physical reconnection model. Mesh/coefficient refinement
and numerical-versus-physical heating must be reported. No theorem of global
monotonicity, positivity or magnetic-shock convergence is claimed.

The nonlinear midpoint iteration evaluates the sensor at its current velocity
and field iterate. Thermal pressure in its signal-speed cap is the entering
gas pressure; no implicit transport closure or extra thermodynamic state is
introduced. A zero coefficient recovers the original centered reference.
