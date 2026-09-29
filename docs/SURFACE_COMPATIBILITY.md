# Thermal compatibility of the sharp resistive surface

The graded study exposes both spatial and time sensitivity. It also motivates
an independent check of the continuum assumptions. The calculation below is a
conditional compatibility test derived for this repository's equations. It is
not a proof of nonexistence of all weak solutions, and it does not identify the
only cause of the measured numerical differences.

## Pressure and traction assumptions

The single-gas pressure equation used here, without conduction, radiation,
viscosity or additional sources, is

```
D_t p + gamma p div(v) = (gamma-1) eta |J|².
```

This follows from the ideal-gas internal-energy equation with physical Joule
heating and mass conservation. The corresponding single-fluid pressure,
induction and Ohm equations are also stated in the
[PPPL M3D-C1 physical model](https://m3dc1.pppl.gov/m3dc1_overview.html).

Across a massless material surface with vacuum outside, no capillary force and
no prescribed surface stress, total magnetic/fluid traction must balance. If
the physical magnetic field is continuous and there is no singular current
sheet, the magnetic tractions cancel and the material-side pressure trace is
zero. This is a conditional statement about the assumed traces. General plasma
interfaces can instead support magnetic jumps and surface currents.
For the traction and magnetic jump conditions, see equations (1.6)–(1.9) in
[Wang and Xin](https://arxiv.org/pdf/2009.11636). Their incompressible,
transversal-field, surface-tension theorem does not validate this compressible
heated, initially tangential-field FRC problem.

## First material-time compatibility condition

For a smoothly transported surface, a smooth pressure trace satisfying p=0
for all time has `D_t p=0` on that surface. If div(v) has a bounded trace, the
pressure equation therefore requires

```
eta |J|² = 0 on the material-vacuum surface.
```

The initialized ideal spherical FRC has radius R, bias B0, and

```
a = -(3 B0/4)(1-s²/R²), s²=r²+z²,
p = (45/8)(B0²/mu0)(r²/R²)(1-s²/R²),
Jtheta = -(15/2) B0 r/(mu0 R²).
```

Its p and B match vacuum without a current sheet, but its volume current has a
nonzero material-side boundary trace away from the axis. Switching on finite
constant eta gives, at initial rest,

```
(D_t p)_surface = (gamma-1)(225/4)(D/R²)(B0²/mu0) sin²(theta),
D = eta/mu0.
```

For the recorded gamma=5/3, R=1 m and D=0.1 m²/s cases, this is 3.75 Pa/s at
the unforced equator and 1.379547904 Pa/s at the initial Gaussian-source
equator, where B0²/mu0=exp(-1). It is strictly positive, whereas persistence
of a smooth zero-pressure trace would require zero. This is a first-time
compatibility defect under the stated smooth, continuous-B assumptions.

The initialized ideal equilibrium is still valid for ideal verification. It is
not, after adding this heating law, an independently known smooth resistive
free-surface solution. The archived moving cases remain useful stress tests;
they cannot serve as proof that a classical heated boundary has been resolved.
Polygonal cells and cell-average pressure do not remove the continuum issue.

## Consequences for development

Retain the original stress inputs and their complete refinement records. Do not
remove edge heating, insert density/pressure floors, or change eta to make the
same result appear converged. Continue the explicit numerical time, mesh,
current and local momentum checks; a continuum compatibility defect does not
excuse an implementation error.

Before using a sharp heated edge for formation experiments, specify which
physical interface assumptions apply. Possibilities require actual justification:
resolved non-smooth expansion, a model for surface currents/stresses, heat
transport, or a density-dependent current response that remains valid as charge
carriers disappear. A mathematically chosen eta taper is not automatically such
a physical closure. Compatible initial data are useful for numerical component
tests, but must not silently replace the original formation problem.

The next physical step is to evaluate a selected species/density/temperature
regime and its current-carrying edge, including collision, Hall and electron
inertia scales. Related interface research discusses electron inertia and
current limitation in dilute wire-array plasma; its geometry, background gas
and numerical density floors differ from this project, so it is a source for
model investigation, not a borrowed validation claim. See Martin's
[Generalized Ohm's Law at the Plasma-Vacuum Interface](https://ecommons.cornell.edu/server/api/core/bitstreams/17448f20-41b4-4457-a472-3ff41ab87529/content),
sections 2.3–2.5 and chapters 6–7.

Spatial prescribed coils and their independent force/work checks remain useful
parallel capabilities. The formation gate must continue to depend on resolved
flow, a defensible thermal/interface model and regime validity, not just those
source capabilities or a conserved energy budget.
