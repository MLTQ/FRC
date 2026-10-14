# Full physical current transport with finite end traces

Retains CurrentTransport's volume advection, changing-density inertia, self-drift,
velocity-gradient stretching and bulk drift stress. Every derivative acts on
the actual full theta and finite axial polynomial/edge current spaces. Interior
radial and axial faces use their two one-sided traces. No face wraps from one
physical end to the other; outer carrier transport flux remains impermeable.

The fluid polynomial is continuous axially, so its interior axial velocity jump
is zero; its discontinuous radial elements retain the original paired stress
tractions. End normal-current consistency follows from the zero toroidal-field
trace in the continuum and still requires a resolved weak-current approximation.
No viscosity, resistive repair or omitted current component is introduced.

Seven full24us coil pulses now exercise this coupled path. Independent limits,
spatial/time comparisons, transport flags and remaining exterior-resolution
limits are recorded in docs/FINITE_CYLINDER_RESULTS.md.
