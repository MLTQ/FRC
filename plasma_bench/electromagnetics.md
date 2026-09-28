# Vacuum electromagnetics
`loop_field` broadcasts coordinates and returns psi, Br, Bz per ampere-turn.
It uses complete elliptic integrals with parameter m and stable small-m/axis
series. Ideal-filament singularities raise ValueError; they are never clipped.
`VacuumField` caches coil geometry bases and superposes prescribed source values
and derivatives. `Field` includes flux, magnetic components, magnitude, dpsi/dt
and Etheta. No plasma solver dependency. Positive current means +Bz on axis;
psi is axial flux divided by 2pi. See PHYSICS.md for equations and limitations.
