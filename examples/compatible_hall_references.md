# Initial states for coupled Hall verification

`wave_state` independently integrates the actual continuum rho, rho v, total
energy, angular momentum and g for a finite-pressure Hall mode. Corner psi is
exactly sampled. It declares compatible magnetic energy but does not adjust E
to make coarse-cell pressure match a target: subcell recovery error remains
measurable. The wave is linear, so finite-amplitude nonlinear residuals remain.
Nonpositive continuum density/pressure amplitudes fail instead of being clipped.

`nonlinear_state` uses the exact manufactured magnetic/carrier snapshot with
zero velocity and prescribed positive ion pressure. Energy is explicitly
constructed under the compatible convention. It is an unbalanced initial
state for real coupled evolution, not an exact time solution or FRC formation.
