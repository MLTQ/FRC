# Finite tube, magnetic vacuum and circular-coil pulses (2026-10-14)

This is physical simulation progress toward frc-vacuum.2, not formation or an
experimental design recommendation. The new solver keeps full finite-mass fluid,
current inertia, Hall force, current transport/stress, two species and Coulomb
exchange. The axial interval now has actual impermeable end caps. The exterior
contains exactly no fluid, conductivity or plasma current; its magnetic energy
extends around the side and both ends to infinity. Prescribed coil filaments
are external sources. No new unit tests, artificial atmosphere, numerical
viscosity, energy repair or residual relaxation was introduced.

## Finite fluid evolution

Five complete nonlinear acoustic trajectories ran. All actual exits0 were
collected (sessions10140,15748). R=.075m,L=1m,n=1e20/m3,Ar+,Te/Ti=5/1eV,
two radial degree5 elements,1600 steps per complete standing-wave period.
The m=3 comparison uses dn/n=1e-5 and the independent linear reflected-wave
solution. Maximum errors relative to perturbation amplitude:

| Axial nodes | Density | Velocity | Pressure |
|---|---:|---:|---:|
|8|0.101785|0.0306781|0.101785|
|12|0.000490840|0.0000800107|0.000490845|
|16|0.0000141367|0.0000148090|0.0000147495|
|16, half amplitude|0.00000706856|0.00000740434|0.00000737492|

Fine-grid remaining departure from linear theory halves with amplitude, as
expected for a nonlinear contribution. Native particle change stays below
1.11e-15 of the total. Energy changes are about3e-5–8e-5 of the tiny acoustic
energy scale, including accumulated floating arithmetic. Stable EOS differences
are recorded separately from absolute-energy subtraction. End pressure force
and discrete zero-velocity constraint reaction are both explicit: maximum
constraint force drops1.81e-5→1.59e-10N under8→16 axial refinement.

The caps represent a sealed tube. They do not model open-end flow, sheath or
particle losses. Full source/equations and states are in finite_tube_* outputs.

## Empty magnetic exterior

The independent field is an off-center axial dipole, alpha=1e-8 Tm3,z0=.1m.
The exterior responds to its prescribed flux trace on the same finite cylinder.
Five matched vacuum patches join a decaying analytic multipole exterior. Outer
trace corrections match the exact analytic trace between nodes, with both
derivatives included in the magnetic energy. Internal vacuum coefficients and
multipole amplitudes minimize that same energy. No material is put in vacuum.

| Radial/axial cells; vacuum cells/modes | Side-field relative RMS | End-field relative RMS | Energy relative error |
|---|---:|---:|---:|
|8/16;8/8|0.283606|0.0117386|-0.118169|
|16/32;12/8|0.0294542|0.000177757|-0.00449996|
|16/64;12/8|0.0134032|0.000143656|+0.0000205500|
|16/64;16/12|0.00446211|0.0000472501|-0.000111701|

Reference exterior energy3.4912538957396843e-7J comes from an independent
160-point Green integral. Energy closeness is nonmonotonic due to cancellation;
the small16/64;12/8 energy error alone does not certify field accuracy. The four
static processes all exited0 (first direct execution and session37631).

## Complete finite-coil evolution

Two radius.12m circular filaments at z=+/-.15m share a prescribed12us sin-squared
current pulse, followed by12us relaxation. Peak on-axis center source B=100uT.
Fluid geometry/initial species are as above; radial fluid elements=field_nr/2,
16 finite axial fluid nodes, field degree5, physical quadrature12. Side and end
flux are independent evolving coefficients. Btheta vanishes on each vacuum
interface through the declared Galerkin space. No periodic coil array is used.

The initial comparison uses fixed eta=1 Ohm m to expose the weak-conductivity
limit Jtheta=-Idot*Aunit/eta. Aunit is independently calculated with a direct
azimuthal Biot–Savart integral, not the production elliptic loop routine.
Expected heat is integral Aunit^2/eta dV * pi^2 Ipeak^2/(2 Tp). This eta is a
mathematical control, not a physical plasma transport recommendation. Finite
induction/inertia, motion and changing species remain in the actual simulation.

Response energy includes fluid/current energy, induced magnetic energy in all
space, and mutual energy I*Phi_ind where Phi_ind=integral_plasma J.Aunit dV.
The prescribed-source response power is I*Phi_ind_dot. Only the divergent
filament self-energy baseline is excluded. Reciprocal linkage is also compared
to vacuum flux evaluated at the actual coil locations. No power error is fitted
or deposited as heat. Accepted Gauss stages supply integrated source work.

Three full fixed-resistance trajectories completed successfully. Actual exits0
were collected for session59968 and both children of batch48323:

| Field grid; steps | Peak current error / reference RMS peak | Final heat / reference | Final energy-work error / reference heat | Coil-linkage discrepancy |
|---|---:|---:|---:|---:|
|8/16;80|0.0334688|0.998896163|1.55334e-5|0.00155491|
|16/32;80|0.00183351|1.000012955|1.55139e-5|0.000322164|
|8/16;160|0.0334575|0.998881282|9.68745e-7|0.00155514|

Spatial refinement improves peak current error18.25-fold. Time halving improves
final energy/work error16.03-fold. The residual current after source turnoff is
not time-resolved uniformly: its final RMS falls5.36e-5→3.38e-6 of the reference
peak when the time step is halved. The Gauss method is not L-stable; the very
short resistive startup/turnoff response is a separate temporal limitation.
Reference heat is7.5258439e-8J. The256→384 azimuthal reference difference is
1.88e-13 relative; production A agrees with the independent integral to2.60e-14
relative. The constant-resistance limit is not an experimental baseline.

Runtime/peakRSS:23.811s/364314624B,67.787s/1177255936B,
44.075s/361234432B in the table order. All states/stages and source closures are
preserved. The initial attempt in finite_coil_r8_z16_t80_eta1 exited1 before
stepping due to a tuple indexing error in vacuum-linkage observation; its frozen
failure remains. The corrected coarse80 case has suffix_v2.

The complete100uT Spitzer pulse finished (actual session98433 exit0,
266.933s,peakRSS1206583296B). It reaches electron Hall parameter0.246 and ion
drift/thermal speed0.792, crossing the existing weak-field/drift review limits.
Peak current20746.6A/m2,peak speed2.01795m/s and final heat55.4676uJ are therefore
conditional model outputs, not validated weak-field transport predictions.
The source self-energy baseline remains excluded; mutual energy is significant
(e.g.-166.58uJ at6us). The current/vacuum linkage discrepancy is0.07011%.

Three complete10uT Spitzer controls stay below those review thresholds. Both
actual child exits0 were collected in batch71975, and actual session74225
completed with exit0. Passing review thresholds is
not a complete transport-validity proof; conduction, viscosity, neutrals and
ionization are still absent. There is no experimental comparison yet.

| Field grid;160 steps | Final heat (J) | Peak current (A/m2) | Peak speed (m/s) | Mean Bz at6us (T) |
|---|---:|---:|---:|---:|
|8/16|5.503183593e-7|1864.40329|0.0198866|2.736674094e-6|
|16/32|5.546933990e-7|2074.21092|0.0201882|2.736811774e-6|
|24/48|5.546878930e-7|2081.10790|0.0201887|2.736810348e-6|

The maximum mean-field difference is1.82557e-5 of the imposed center peak, but
peak current changes10.12% relative to the finer grid and heat0.789%. Thus
average-field agreement does not establish accurate current layers. Fine-grid
Hall/drift ratios peak at0.024596/0.079201. Native particle change stays below
4.74e-17 of total particles; maximum energy/work discrepancy is9.50e-8 of the
stated field-energy scale (not a pointwise field error). These are full moving
plasma runs, with two species evolving to mean Te/Ti4.413688/1.586313eV at24us.

The24/48 trajectory also completed all160 steps and24us. Relative to16/32,
peak current changes0.3314%, final heat0.000993%, peak speed0.00267%, and the
maximum mean-field difference is4.34234e-7 of the center source peak. This is
observed convergence of those quantities, not a uniform error certificate.
Its maximum energy/work discrepancy is9.47436e-8 of the stated magnetic scale;
maximum particle change2.59e-17. Hall/drift ratios peak at0.024622/0.079468.

The independent vacuum-flux/current-reciprocity linkage discrepancy changes
0.07014%→0.07518% on16/32→24/48: fixed12-cell/8-mode exterior resolution is not
jointly converged by refining only plasma cells. That remaining exterior
approximation is explicit; energy closure does not remove it. The static dipole
study above separately demonstrates exterior-resolution improvement. No extra
coupled exterior-resolution campaign is claimed here.

At source peak the plasma-volume mean Bz is2.73681uT, versus9.50503uT for the
coil-only field averaged over the same volume (28.79%). Finest peak current is
2081.11A/m2; final thermal gain is0.554688uJ. All these outputs describe the
conditional fully ionized, sealed10uT model, not a proposed formation shot.
Runtime/peakRSS for8/16,16/32,24/48:45.999s/351780864B,
137.981s/1121828864B,354.970s/2753855488B. Every run remained within its4GiB
self guard. No scientific worker remains live.

The saved-field figure was visually inspected; its observer69616 exited0.
The figure is outputs/finite_coil_comparison/finite_coil_response.png
(and PDF from the same figure). It compares16/32 and24/48 histories and shows
actual Jtheta, Bz and poloidal flux contours at6us. It is a new pulse figure;
the older interactive3D viewer still shows its explicitly labeled old contact.

## Remaining boundary

The bounded sealed finite-coil component is implemented and its measured
comparisons are complete. Continue frc-vacuum.3: moving free surfaces and
physically supported end/transport
conditions, ionization/model validity and actual formation. No experimental
baseline, field reversal, closed-flux formation,3D stability or confinement is
established by these finite-tube comparisons. The full FRC goal remains open.
