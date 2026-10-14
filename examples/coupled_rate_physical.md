# Physical full-rate power and cross-grid response reconstruction

`response` independently reconstructs the centered physical RHS at the actual
signed states, including both EOS pressure derivatives, Jdot and algebraic E.
It returns the same field fingerprints used by the frozen tangent audit.

`power` evaluates physical instantaneous fluid, current and magnetic energy
derivatives with actual reconstructed n,v,J,B and saved full returned rates:
fluid power is M(nt |v|²/2 + n v·vt) + sum pt/(gamma-1);
current power is mred/e² (J·Jt/n - |J|² nt/(2n²));
magnetic power is B·Bt/mu0. Actual input current coefficients are supplied from
the unchanged weak Ampere map, not inferred by floating amplitude scaling.

Particle, axial momentum, physical angular momentum and both entropy-density
integral rates are separate. Total power is normalized by the sum of absolute
physical powers; there is no residual subtraction or normalization floor.
These instantaneous global measurements do not close local or time-integrated
enthalpy/Poynting/vacuum budgets. All units and finite-amplitude effects remain.
No evolution, solver repair or tests.
