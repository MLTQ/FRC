# Numerical compression treatment and the long-array pulse

This is an optional solver change for strong compressive flow, not a physical
transport coefficient or an established FRC formation model. The original
nonviscous trajectories remain unchanged and preserved. The option defaults off.

## What the actual strong pulse exposed

The independent nine-loop20mT pulse failed at1.064us with1.75ns RK4 steps;
its last accepted state was already numerically explosive. The split-current
method also failed after613steps. A complete restart with0.875ns steps instead
reaches1.101625us before a stage folds a material cell. A short0.4375ns continuation
from the earlier step607 reaches a similarly skewed cell near an end. It does
not establish full time convergence because it shares the earlier coarse prefix.

Saved original stages show field/current growth before the entropy explosion;
Coulomb temperature exchange is much slower on that step. The restricted
fixed-fluid EM eigenvalues do not demonstrate a simple pure electron-gyro
stability limit. Full-system stability and spatial accuracy remain unresolved.
Read outputs/material_independent_pulse_failure_analysis and the complete failed
outputs/material_tube_independent_20mt_r4_z16_t4000; no states were repaired.

## Explicitly numerical pressure and heat

For theta=div(v), c²=sum(gamma_s*p_s)/rho and minimum physical cell altitude h:

    u = h max(-theta,0)
    q = rho (0.5 c u + 2 u²)
    force = integral q div(delta x) dV
    heating = -q theta >= 0

The same cylindrical geometry/quadrature pairs the force and heat. The force
changes centre-of-mass momentum; the declared heat partition puts the loss in
ions through the actual dual entropy representation. Numerical heat is reported
separately from eta J². It does not alter physical resistance or add relative
current drag. Exterior matter/conductivity stays exactly zero. In smooth flow
q decreases as O(h)+O(h²), and q=0 in expansion. This scalar treatment cannot
resolve all shear/hourglass modes or guarantee an untangled mesh. Energy after
finite RK4 steps remains an observation, never a fitted correction.

The method uses the artificial-pressure idea of von Neumann and Richtmyer,
[1950 paper](https://doi.org/10.1063/1.1699639), and conjugate force/work reasoning
discussed in [Caramana et al.1998](https://cnls.lanl.gov/~shashkov/papers/comp.pdf).
It does not implement that paper's complete staggered algorithm or claim its
exact fully discrete conservation; the formula above is the implemented choice.

## Complete independent wall-shock baseline

The same material fluid, mass and dual heating evolve a uniform gamma5/3 gas
with Ar ion plus electron mass, n=1e18/m³ and10km/s incoming flow. Initial total
pressure is0.01 rho u0². A fixed end wall stops the downstream fluid; rigid radial
constraints give a one-dimensional shock. The moving upstream boundary's reaction
work is integrated using the actual RK stage weights. No electromagnetic or
species exchange terms belong to this Euler fixture.

Mass/momentum/energy jumps give shock speed3457.183m/s toward the inflow,
compression3.89252846 and p2/p0=135.5718293. At20us the exact front isz=0.43085634m.
These values are independent of solver output. All three complete trajectories
finish with the following volume L1 errors:

| Axial cells | dt, ns | Density | Velocity | Pressure | Energy/work residual / E0 |
|---:|---:|---:|---:|---:|---:|
| 50 | 50 | 3.1206% | 1.6184% | 10.4266% | 0.0045% |
| 100 | 25 | 1.9824% | 0.9162% | 6.1963% | 0.1049% |
| 200 | 12.5 | 0.9623% | 0.4845% | 2.9564% | 0.0523% |

The inspected PNG (also saved as PDF) shows narrowing shock width and correct downstream
pressure, with a localized density/wall-heating defect still visible. The first
50-cell mesh uses its radial altitude as h; on finer meshes h is axial. The
energy residual is not monotone from50 to100, and the given timestep refinement
alone does not remove heating-integration error. A complete100-cell run at four
times as many steps (3200steps,6.25ns) reduces the relative energy/work residual
from1.04936e-3 to3.08385e-6 (340-fold), with density/velocity/pressure L1 errors
1.8912%/0.9074%/5.7313%. It retains the correct shock profile. The inspected
PNG and generated PDF are in outputs/material_wall_shock_viscosity_t4, source
be052c85250c68fe201c8832b6e4d2e6ea3a5940b101b28a3a51dbb93745f030.
This isolates a substantial timestep contribution; the remaining wall-density
defect and shock width are not removed by that time refinement.

Evidence: outputs/material_wall_shock_viscosity_b; frozen source
112a2fde261fde404b4de8c1cd0755a37f330af0cb60acbb88bc1938df26365b.
The earlier single-radial-cell setup failed the mesh contract before evolution
and remains in outputs/material_wall_shock_viscosity. No unit tests were used.

## Physical pulse in progress

LIVE13616: outputs/material_tube_independent_20mt_viscous_r4_z16_t4000.
Same complete initial plasma, independent bias/long-array waveforms and0.875ns
step, now with C1=0.5,C2=2 numerical bulk pressure. Frozen source
05555dc9be49c3a1af9c31467fc713f4d5e395442328e9107eea19f68d0b3c75.
This run has not yet established whether the end fold is resolved. Collect its
actual terminal state; preserve failures and do not edit its source while live.
A completed damped trajectory would still need physical comparison and spatial
assessment. It would not establish ionization, transport validity, confinement,
three-dimensional stability or a laboratory design.
