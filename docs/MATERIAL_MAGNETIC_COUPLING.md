# Moving magnetic force and winding (2026-10-14)

This extends the material fluid from 05b4abd4 with actual magnetic forces and
angular torque. It is a necessary coupling component, **not the completed
finite-mass moving plasma solver**. The full finite-coil/current/species solver
still has a fixed footprint; this moving component currently uses ideal MHD.
The full original goal and frc-vacuum.3 remain open.

## Equations and boundary

The evolved canonical variables are material positions, meridional momentum,
angular momentum, material rotation angle, and two specific entropies. Existing
consistent physical mi+me inertia, adiabatic pressure forces and centrifugal
force are retained. The full cylindrical deformation gradient F contains both
meridional deformation and r*grad(theta). Its determinant is the actual volume
ratio. The magnetic field follows the ideal Cauchy formula

    B = F B0 / det(F),
    Wb = integral dV0 |F B0|² / (2 mu0 det(F)).

This is the frozen-flux ideal limit of Faraday/Ohm and MHD momentum, as described
by the [PPPL physical model](https://m3dc1.pppl.gov/m3dc1_overview.html) and
[Fitzpatrick's flux-freezing derivation](https://farside.ph.utexas.edu/teaching/plasma/lectures/node106.html).
The implementation differentiates Wb analytically with respect to both position
and material angle. These forces and torques enter the same RK4 trajectory as
pressure and inertia. Positions and angle are unconstrained apart from the
regular symmetry axis. Full poloidal/toroidal winding is retained. There is no
prescribed shape, force fit, numerical resistivity or energy repair.

The reference flux is psi=(Bp/2) r0² (1-s²)^4, with
Btheta0=Bt*(r0/R)*(1-s²)^3 on the actual reference material mesh. The magnetic
field vanishes at its boundary. Its ideal pullback therefore permits exactly
zero exterior B and matter: no exterior atmosphere or conductivity, no wall,
no coil and no imposed magnetic boundary stress. The initialized closed flux
is part of the fixture; it is **not newly formed flux or reconnection**.
This zero-field exterior is a restricted boundary contract, not the general
moving plasma/vacuum/coil boundary required for the laboratory sandbox.

## Current reconstruction issue found and corrected

Q1 geometry and material field gradients permit tangential B jumps between
cells. A current reconstructed only by differentiating within each cell misses
these distributional face contributions. The first complete r4 run exposed
this: its initial normalized current error was 0.4103, despite much smaller B
and force errors. That trajectory and its original source remain preserved.

The current display now solves the physical weak Ampere relation

    integral G.J dV = integral curl(G).B dV / mu0,

using regular vectors (r*N,0,0), (0,r*N,0), (0,0,N). The outer boundary term
vanishes because B=0 there. This projection includes internal face contributions;
both it and the within-cell curl remain recorded. Dynamics are unchanged by
this diagnostic correction. The projected current is not an adopted evolving
carrier state, a pointwise charge-conservation proof, or a guarantee of bounded
J/n at the vacuum edge. Its inferred current energy and drift diagnostics are
representation-dependent. They cannot justify dropping finite-mass terms.

## Independent physical comparisons

The spherical continuum fields are differentiated independently of the material
operator. With h=1-(r²+z²)/R²,

    Br = 4 Bp r z h³/R²,
    Bz = Bp h⁴ - 4 Bp r² h³/R²,
    Btheta = Bt r h³/R,
    Jr = 6 Bt r z h²/(mu0 R³),
    Jtheta = 4 Bp r h² [5-11(r²+z²)/R²]/(mu0 R²),
    Jz = Bt [2h³-6r²h²/R²]/(mu0 R).

These give independent initial fields/current and JxB acceleration. Pressure
adds C*(r,z)/R², C=5e*(Te0+Ti0)/(mi+me). The actual polygonal reference differs
from the exact sphere, and that initial geometry error is measured. A separate
80x80 spherical quadrature supplies the continuum magnetic energy.

The all-space MHD scalar virial identity supplies a complete trajectory
comparison independent of the force assembly:

    I = integral rho |x|² dV,
    I'' = 4 K + 6(gamma-1) U + 2 Wb.

For gamma=5/3 and conserved E, this integrates to

    I(t)=I(0)+I'(0)t+2 E(0)t²-2 integral_0^t (t-s) Wb(s) ds.

The observer uses a sampled trapezoidal integral and includes rotational kinetic
energy. Nothing is fed back to forces. Conservation and virial agreement alone
do not establish local field accuracy; the spatial comparisons remain separate.

## Scope still missing

The run is a formal ideal-MHD fixture. Finite current inertia, Hall response,
current advection/stress, species relative entropy transport, Joule heating,
Coulomb exchange, neutrals/ionization and finite coil work are not implemented
in this moving component. It must not replace the existing fuller fixed-footprint
solver. The cold dilute edge remains outside a demonstrated physical closure.
The field deformation/energy derivatives built here are intended for that full
coupling. Do not launch more ideal campaigns as a substitute for it.

Next implement independent evolving magnetic/current coordinates on the moving
geometry, retaining the full current-energy geometry work and both species
transport. The weak-current result shows why cellwise curl and pointwise J/n
cannot simply be bolted onto this Q1 field. Couple the existing magnetic-only
vacuum and finite coils with their true shared flux traces and source work.
No ideal closed-flux result establishes FRC formation, stability or confinement.

## Completed trajectories and measured limits

Each run completes800 steps through7.1217692648 microseconds: R=35mm,
Ar+ mass, nc=1e24/m3, central Te=Ti=1eV, Bp=0.3T, Bt=0.15T. These are
mathematical ideal-MHD inputs, not selections for the user's laboratory fill.
The original r4 run (actual session95584) and three weak-current-observer runs
(batch18833) exited0; every accepted state and complete source are saved.
Figure observer99255 also exited0. No scientific worker remains live.
No unit tests were added or run.

Initial RMS errors are normalized by Bp, Bp/(mu0 R), and C/R respectively:

|Grid|B error|J error (weak)|Acceleration error|Initial Wb relative error|Final Wb [J]|Final peak vtheta [m/s]|
|---|---:|---:|---:|---:|---:|---:|
|4/8|0.02266014|0.13819156|0.016952233|-0.027283341|0.15093541|4.4305167|
|8/16|0.010648532|0.045418992|0.004331837|-0.0069295401|0.15429227|4.4691947|
|16/32|0.0052379098|0.01163322|0.0010637349|-0.0017391697|0.15515263|4.5396801|

Acceleration improves about fourfold per mesh doubling, and weak-current error
improves3.04x then3.90x. The field error only improves2.13x then2.03x; do not
claim second-order local B accuracy. Native energy change is below2.89e-14;
virial inertia error stays below6.38e-10 of initial inertia. These do not erase
the initial magnetic-energy error, still0.174% on the finest polygonal mesh.
Within-cell current error instead approaches0.3613 of Bp/(mu0 R), illustrating
why it is an incomplete current diagnostic.

Final common-material-label RMS differences (not Eulerian errors) from4→8 and
8→16 cells are respectively: position0.0063451/0.0014972 of R;
velocity0.0037171/0.00090620 of sqrt(C); B0.0123093/0.0061582 of Bp;
J0.0616262/0.0147149 of Bp/(mu0 R); density0.0034153/0.00112090 of nc;
angle7.33883e-5/1.74990e-5rad. Finest final peak B=0.134018T,
peak projected |J|=9.30935MA/m2 and peak toroidal speed=4.53968m/s.
The full toroidal motion was initially zero and is generated by magnetic torque.
These scalar/common-label comparisons do not establish a high-accuracy
resistive formation trajectory or convergence of unresolved face currents.

The inferred edge drift/ion-thermal ratio grows under refinement (final
6.39,48.2,78.3); ion skin depth also exceeds the body radius in dilute samples.
The continuous current projection has no enforced carrier-weighted edge trace,
so these large ratios mix numerical reconstruction and physical model limits.
They are not a validated transport calculation. Inferred current energy is
small globally but cannot certify pointwise ideal-MHD applicability.

The inspected scientific figure is
[material_magnetic_expansion.png](../outputs/material_magnetic_comparison/material_magnetic_expansion.png);
its PDF and common-label comparison JSON are beside it. Initial/final panels
share physical color scales. The previous interactive3D viewer is unchanged.

Original source9168b7d270eb3c8ce955e0b3903576b860264b8c95754b57f64d0a44dd3d5064;
three refinement sources1b894fc1384b4bad93cc5d9f7018f90d3f5231617749a7fe010791263482a38b.
Frozen sources remain authoritative; the diagnostic correction did not modify
the magnetic-force dynamics. Scientific and figure workers are terminal.
