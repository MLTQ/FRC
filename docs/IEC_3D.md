# Positive-cage 3D kinetic model

The first electrostatic kinetic path is built and runs locally through the
[experiment app](EXPERIMENT_LAB.md). It is a collisionless development model for
frc-cac.2, not the user's historical air discharge and not a fusion/POPS prediction.
The goal remains to discover whether a self-consistent electron population can
produce a useful ion well and a driven response despite losses.
The user's adaptive-pulse variant is **HOPS — Harmonically Oscillating Plasma
Sphere**. Its first actual +3 kV / 0–5 µs experiment and matched controls are
documented in [HOPS.md](HOPS.md).

## Equations, geometry and sources

Electrons and singly charged ions obey dx/dt=v, m dv/dt=q E with actual masses;
the default ion mass is 2 proton masses (a hydrogenic reference, not molecular
chemistry). At every step, solve -epsilon0 Laplacian(phi)=rho, E=-grad(phi).
There is no imposed central well, symmetry sector, neutralizing background,
density floor, temperature rescale, particle reflection or energy repair.

The gas domain is a grounded spherical chamber of radius 20 cm, excluding the
full 15 cm truncated-icosahedron frame: 90 circular 2 mm wires and 60 rounded
junctions. The sphere is an explicit reference assumption; actual chamber shape,
feedthrough and supports are not represented. Steel is equipotential and absorbs
particles at their first mesh-face crossing. No secondary electron emission.

A specified radius-5-cm seed contains independently sampled 3D electron/ion
positions and Maxwellian velocities. Globally equal initial species counts do
not enforce subsequent neutrality. The optional source emits physical electrons
from uniformly sampled grounded-wall area, monoenergetically along the inward
surface normal. Its current is a Poisson mean; actual charge/energy are counted.
This is an idealized source experiment, not air ionization or a specified emitter
design. Bias V=V0[1-exp(-t/rise)][1+depth sin(2 pi f t)], with rise=0 selecting DC.
No resonance is selected by the force law.

Gmsh builds a conforming tetrahedral gas mesh, fine near wires. Linear FEM uses
Kij=epsilon0 integral(grad Ni dot grad Nj); qi=sum_p qp Ni(xp). A float64 sparse
Dirichlet solve is prefactored. E is the gradient of the same interpolated
potential used for particle shapes. Kick-drift-kick integration traces each
straight drift through tetrahedral faces to its first conductor hit. Emission
precedes the full step's first kick. No explicit numerical diffusion/filter is
added; planar metal, finite particle shape/noise/self-force and piecewise constant
E remain truncation errors. This path currently has no GPU acceleration.

## Conservation

Electric field energy is phi^T K phi/2. Metal charge is (K phi-q) on each
Dirichlet boundary. Supply work integrates V dQ_supply with
dQ_supply=dQ_electrode-q_collected; grounded-wall emission contributes injected
kinetic energy. The uncorrected ledger is K+field+absorbed-injected-supplywork
minus its initial value. Particle kinetic throughput is also reported so the
much larger electrode charging energy cannot hide particle-scale errors.
Charge lost to each electrode and actual emitted charge close the particle
charge ledger. Finite dt boundary kicks and midpoint supply work have error.

## Independent measured evidence

`outputs/iec_fem_analytic_v2/` compares uniform volume charge in a grounded sphere
with Gauss law. At 40/28/20/14/10 mm targets, potential L2 errors fall
1.817/0.873/0.434/0.215/0.110%; volume-weighted field errors fall
6.106/4.082/2.856/1.963/1.389%. The independent concentric-capacitor field errors
at 40/28/20 mm fall 21.77/15.95/11.86%, and capacitance/field-energy errors
5.82/3.12/1.73%. These are convergence evidence, not near-wire certification.

The same sphere's exact electron oscillator exposes coarse force error: a small
3.44 cm-amplitude orbit over one analytic period still has 12.6% endpoint error
and 5.05% maximum analytic-energy error at 10 mm / 200 steps. Spatial phase
errors are not monotonically converged at every level; reducing dt alone does
not fix piecewise constant E. Small-well orbit accuracy remains a real gap.

`outputs/iec_cage_fem_refinement_v1/` refines the actual cage:

| Bulk / wire target | Nodes | Capacitance | Centre phi / V |
|---|---:|---:|---:|
| 24 / 1.5 mm | 28,917 | 11.5707 pF | 0.884371 |
| 18 / 1 mm | 68,353 | 11.3948 pF | 0.878958 |
| 14 / 0.75 mm | 127,802 | 11.3041 pF | 0.876060 |

The separate boundary method gives 11.090 pF and 0.866075, with up to 1.4% local
surface-potential error itself. Cross-method differences shrink; neither is an
exact local-gradient or wire-loss reference.

The app-launched `outputs/lab_runs/iec_20261021_233257_b1cb/` and
`outputs/iec_3d_static_dt10_20ns/` are paired actual 20 ns, 100 V DC, source-off
PIC controls with 600 particles/species on the 18/1 mm mesh. At nominal 20→10 ps,
endpoint position difference is 0.00223%, velocity 0.0466%, electron kinetic
energy 0.0111%. Absolute ledger residual falls 1.29e-15→1.82e-16 J, or
0.00924→0.00130% of particle kinetic energy. Charge closes to roundoff. These
controls support timestep consistency, not particle/mesh convergence.

## Current physical scope and next action

Actual source-driven 0–50 ns runs use +100 V with a 5 ns rise, 100 µA ideal wall
source at 2 eV, n0=1e11 m^-3, Te=1 eV, Ti=0.05 eV, 600 macros/species and the
18/1 mm mesh. Initial Debye length is 23.5 mm. These are declared exploratory
inputs, not an optimum. The final app-launched run
`outputs/lab_runs/iec_20261021_234301_0adb/` completed 0–50 ns with 101 actual
frames and full 3D state archives in 117 s. It ends with 893 electron macros,
600 ions and 55 electron cage hits. Maximum absolute energy-ledger residual is
1.42e-14 J (final 0.0106% of particle kinetic/injected/absorbed throughput), and
maximum charge closure error is 1.22e-26 C. Centre potential is +88.23 V at
+100 V cage bias; this scalar is not a trapping or POPS result. Source hashes,
all frame times/populations, final raw fields/energies and charge deposition
were checked in `outputs/lab_app_delivery_v1/report.json`.

No neutral/Coulomb collisions, ionization, radiation or fusion reactions exist.
Initial/evolving local resolution, source distribution, particle noise and
near-wire losses need convergence; 50 ns cannot establish ion oscillation,
resonance, compression gain or long-term trapping. Next: improve and measure
small-well force accuracy and particle/mesh sensitivity, then run electron-well
formation, free ion response and driven/DC comparisons over actual ion periods.
Only after physical losses/chemistry are represented should pressure/source/drive
scans be interpreted as favorable experimental parameters. FRC jobs continue.

The HOPS campaign completes actual source-off +3 kV / 20 ns comparisons:
10→5 ps lowers maximum particle-throughput energy residual 0.378→0.113%, with
charge error below 1.87e-27 C. Its actual 0–5 µs feedback run stopped at
2.52 µs for an unknown reason. On Oct 25, source-stream recovery and an actual
2.51→2.52 µs solver replay matched all particle/field arrays bitwise, physical
diagnostics exactly, and source RNG/controller state. After recovery, the user
prioritized pulsing/relaxation; v1 is retained at a full checkpoint at 2.628245 µs
with no gate switches, and its DC/fixed tail is stopped. The next probe is one
0.2 µs +3 kV command followed by source-free relaxation through 5 µs, with
20 ns electrode response and unchanged 5 ps steps. See [HOPS.md](HOPS.md),
`outputs/hops_pulse_release_20261025_v1/` and the original restart evidence in
`outputs/hops_resume_20261025_v1/`. Spatial/phase and
physical collision losses remain unresolved; this is a development experiment.
