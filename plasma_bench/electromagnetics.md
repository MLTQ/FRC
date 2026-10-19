# Vacuum electromagnetics
`loop_field` broadcasts coordinates and returns psi, Br, Bz per ampere-turn.
It uses complete elliptic integrals with parameter m and stable small-m/axis
series. Ideal-filament singularities raise ValueError; they are never clipped.
`VacuumField` caches coil geometry bases and superposes prescribed source values
and derivatives. `Field` includes flux, magnetic components, magnitude, dpsi/dt
and Etheta. No plasma solver dependency. Positive current means +Bz on axis;
psi is axial flux divided by 2pi. See PHYSICS.md for equations and limitations.

## Shared-waveform grouping and time cache (2026-10-19)
`VacuumField` sums coils that share one waveform object once at construction.
A coil's current is turns × polarity × waveform(t), so each group contributes
waveform(t) times a fixed weighted basis. `at(t)` therefore evaluates each
distinct waveform once instead of looping over every filament, and keeps the
last 8 times in a cache.

**Effect:** in the 112-filament NUCTE-III setup, the step cost fell from
75 ms to 14.7 ms. Fields agree with the per-coil sum to about 1e-15 relative
(summation order only).

Returned arrays are read-only, so in-place modification by a caller raises
instead of corrupting the cache.
