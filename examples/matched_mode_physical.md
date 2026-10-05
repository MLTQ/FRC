# Independent matched-mode physical reconstruction

Root-polynomial and independent Fourier quadrature reconstruct physical density,
velocity, species specific entropy/pressure/temperature, current and magnetic
fields. `signal` extracts the centered physical difference of actual opposite
perturbations, including nonlinear quotient and exponential equations of state.
No continuum mode is projected onto production eigenvectors.

`changes` independently integrates the original bulk kinetic, species internal,
physical current-inertia, magnetic increment and harmonic bias cross energies.
Stable conservative differences/log1p/expm1 avoid subtracting large thermal
backgrounds; energy/state is never repaired. Particle and both log-entropy
changes, point/sample density and temperature extrema stay separate. Spatial
quadrature convergence is required; this is distinct from grid conservation.

Background bias has no current. Current is reconstructed from original weak
Ampere coefficients, with physical masses; this is an energy audit, not a claim
that all strong current/interface balances vanish. The caller supplies explicit
positive epsilon, identical meshes and actual saved signed states.
