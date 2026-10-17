# Density interfaces and mechanical resolution

## All-heavy captured stage completes; full coupled comparison underway

The original all-heavy captured 5 ns stage is terminal successful (81727 EXIT0),
350 steps, all 227,649 density cells retained, Erel 3.29120e-12. Ar5–Ar9 endpoint
peaks fall from 4.21/5.63/7.87/14.18/370.16 eV to .653/.698/.778/.972/.468 eV.
Ar8 volume RMS excess above reversible compression falls .722809→.0234416 eV;
Ar9 falls 29.7306→.0157336 eV. Full-support and interface plots now include
original, Ar9-only and all-heavy endpoints and have been inspected.

Lower-charge edge cooling increases, and their volume RMS excess can increase:
Ar0 .00000355→.00078439 eV. This is not general thermal convergence. An adiabatic
reference is not an exact shock solution, and no population is omitted because
it is small. Independent 80-digit widths reproduce the new Ar8/Ar9 reversible
reference within 5.55e-17 eV. Full data and the all-species table are in
outputs/reaction_density_jump_heavy_transport/README.md. PID6975 is finished.

Eliminating duplicate within-call width calculations preserves every saved
array in a complete matched 20 ns all-heavy/10 ns-source run: 1,024 steps,
Erel 2.63129835904e-11. This is a runtime simplification, not a physical change.
The new full 20 ns all-heavy/5 ns-source trajectory uses that evaluated source
(7035e61481459796986bded32038e00a468d8133563ea45d5c17c52fbfdb3684), frozen under
/private/tmp/frc-all-heavy-coupled-5ns-widths-20261017, session30896.
It is now terminal failed (30896 EXIT1): completed coupled time 10 ns, followed
by 73 accepted internal steps in the third interval before the acoustic
generalized eigensolver reported '2 eigenvectors failed to converge'. Source
moments, exact inherited knots, accepted steps and failure are preserved. This
is not a full coupled success. That state now has an exact73-step replay and independent100-digit spectrum;
see ACOUSTIC_MASS_SCALING. The corrected run completes third-interval transport
(2058steps/full5ns) but rejects thermal energy in the following source half.
Its reaction_source_gather_failure reconstruction is live as2255. Extra
enrichment also prevents a pure temporal-
order claim. All spatial/current/physical-closure requirements remain.

## Complete coupled20ns candidate

The original Ar9-only density-interface trajectory now reports terminal success
at20ns and its process has exited. It takes972 transport steps, retains all
2,414,399 density cells and has Erel2.79334e-11. Final Ar9 maximum falls from
370.173188550eV in the matched scaled-heavy5ns baseline to0.637746842eV.
Ar8 remains4.432301209eV, essentially unchanged. The last Ar9 source return
changes T by at most0.009365589eV. These are complete coupled endpoints, distinct
from the isolated5ns stage comparison below; the two refinements start at
different points in the trajectory and must not be conflated.

Raw endpoint/source/actual-step evidence is in
outputs/reactive_material_argon_density_jumps_5ns. Its completed precise-union
comparison includes455,321 intervals, sampled in8192-interval batches to bound
memory. E/J L2 differences relative to the original trajectory are
1.0999492e-7/1.6886118e-4. Ar9 common-volume temperature and velocity RMS changes
are29.675866eV/6995.04289m/s; Ar8 temperature RMS changes0.003239458eV and
electron T RMS7.40178e-5eV. Front differences are1.87180e-12m for Ar9 and
3.65177e-11m for electrons. Pointwise maxima at displaced discontinuities remain
recorded. This quantifies the local correction without inferring global accuracy.

The evaluated source-return controls are now integrated as explicit options,
disabled by default. Two fresh full20ns all-heavy10ns-source trajectories from
main now complete, with contrasts10000 and100. Both take1024steps and retain
21,385density cells with Erel2.6313e-11. Against the original10ns-source control,
contrast10000 changes E/J by8.46607e-7/3.24718e-4 relativeL2 and Ar8/Ar9 maxima
become0.202018/0.201416eV. Lower-charge edge cells cool to0.0078–0.0145eV;
the neutral center stays0.1eV. Neutral/Ar1 front differences are42.24/41.74um.

Contrast100 actually adds more velocity nodes and changes E/J from10000 by
5.04981e-15/3.05631e-6. Ar8/Ar9 common-volume T RMS changes2.26503e-4/3.21265e-5eV;
their velocity RMS changes31.6235/6.86396m/s. Strong edge cooling persists.
These are threshold-sensitivity measurements, not proof of thermal/space/time
convergence. READMEs and complete raw evidence are in the all_density_jumps_10ns
and all_density_jumps_c100_10ns output directories. PID6975/session81727 is now terminal complete
on the preserved all-heavy captured5ns stage (see the later evidence above). Physical equations, masses,
exact vacuum and represented populations remain unchanged.

## Completed actual-stage spatial comparison

`outputs/reaction_density_jump_transport/comparison.json` and the inspected
`interface_stage.png`/`interface_detail.png` compare the full original 5 ns
endpoints, retaining all 21275 Ar8 and 15033 Ar9 material cells. Independent
80-digit Decimal differences of saved positions reconstruct the reversible
reference T0*(h0/h1)^(2/3); the maximum difference from the saved references is
2.22e-16 eV. This reference is isentropic compression, not an exact shock
temperature: physical shock entropy can legitimately exceed it.

| Ar9 endpoint quantity | Original 38 velocity elements | 42 elements with density interfaces |
|---|---:|---:|
| Maximum T (eV) | 370.160809454 | 0.468429806 |
| Maximum excess over adiabatic (eV) | 369.748653523 | 0.244441032 |
| Volume RMS excess (eV) | 29.730557863 | 0.015733633 |
| Particle RMS excess (eV) | 4.47845e-9 | 1.89336e-9 |
| log10 density at the hottest cell (m^-3) | -293.780775 | -293.996959 |
| Right support edge (mm) | 8.674128406 | 8.673367130 |

The original hottest cell is at 8.2034 mm; the resolved hottest cell is at
7.9930 mm. Particle-weighted and global energy measures conceal substantial
local errors in these extremely dilute retained populations. No population is
discarded on that basis. Ar8 remains essentially unchanged: maximum 14.176049
eV at -9.9578 mm, versus 0.910445 eV reversible reference there, with volume
RMS excess 0.722809 eV. This supports the existing all-heavy interface trial.
These are final-stage extrema, not maximum-over-time measurements or a claim
that a laboratory discharge produces these high charge states. No threshold,
spatial or source-time convergence is established by these two endpoints.

Latest execution/result note (2026-10-17, 13:49 UTC): the full Ar9-only stage
16298 is EXIT0 at 5 ns after 221 steps, Erel 5.1023017e-12. Its summary reports
final Ar9 maximum 0.468429806 eV and maximum/RMS excess over adiabatic
0.244441032/0.015733633 eV. The matched original control ended at 370.160809 eV;
Ar8 is unchanged near 14.176049 eV. This does not establish general convergence
or fix the other heavy species. Raw outputs are preserved; complete coupled
and all-heavy comparisons remain required. After user-authorized tool shutdown
and uv cache cleanup, coupled PID 2415 resumed; heavy PID 6975 remains suspended.
Wall time includes the long storage pause. See the latest THREAD_HANDOFF.

Earlier checkpoint:

Current execution note (2026-10-17, 13:28 UTC): cache cleanup recovered about
6.3 GB. Ar9-only16298/PID1882 was identity-checked, resumed and confirmed running;
coupledPID2415 and all-heavy81727/PID6975 remain suspended for serial scheduling.
Source/data are unchanged. See THREAD_HANDOFF and
outputs/storage_cache_cleanup/cleanup.json. Collect the first full endpoint,
then verify identity and available export/swap space before resuming each
remaining process. Earlier storage failures remain preserved.

The active lab-facing goal remains unfinished. This experiment addresses the
local thermal-transport defect in `frc-transport.5.5.5`. It is a concrete
neutral Euler comparison followed by actual argon-stage and coupled-run trials;
it is not a formation experiment or a new physical heating model.

## Independent physical problem

Both runs use neutral argon with physical mass39.948u, initially zero velocity,
T=0.2eV, n=1e14m^-3 to the left of zero and1e6m^-3 to its right. The support
is[-10,10]mm, exterior vacuum is exactly empty and the remote boundaries are
at±0.2m. Both request the original1.1510438318865044us duration, before central
waves interact with the exterior vacuum fans. The solver retains P0 density/
thermal energy, P1 velocity, consistent physical inertia, gamma5/3, minmod
numerical viscosity(.5,2) and half of its actual acoustic bound. Numerical
viscosity is explicitly numerical entropy production, not collisional heating.

The same65 fine density cells represent both initial continuum states. Shifting
the interior nodes of a64-element uniform velocity mesh by0.99 element widths
puts the density jump inside a velocity element, with the dense material in
only1% of its width. The resolved run adds that one existing fine density face
to the velocity basis. It keeps all old velocity nodes, fine cells and amounts.

The exact central solution follows the independent Rankine–Hugoniot shock and
isentropic rarefaction relations in `examples/material_density_step.md`. Its
maximum temperature is1.2042418371603956eV. Independent exterior free edges
are x=±(a+2*c0*t/(gamma-1)), c0=sqrt(gamma*e*T0/m), giving edge speed
2691.8059861263746m/s. The exterior fans cool and cannot raise the maximum.

| Measured quantity | Jump inside element | Node at jump |
|---|---:|---:|
| Actual velocity elements |64|65|
| Fine density cells |65|65|
| Requested duration(us) |1.151043832|1.151043832|
| Accepted duration(us) |0.637895323|1.151043832|
| Accepted steps |4655|610|
| Full-history maximum T(eV) |3040.625629|0.925315439|
| Full-history maximum speed(m/s) |334121.929|1873.571|
| Maximum front error(m) |0.187593522|0.002111746|
| Maximum relative energy error |1.074e-11|3.635e-11|

The buried case fails when its right material edge approaches0.2m. Its failed
prefix is not a full endpoint. Its temperature peak occurs at15.05585655ns;
looking only at its final central temperature(~0.205eV) hides the runaway.
The resolved run finishes, but its final central peak0.909386556eV is below
the exact1.204241837eV, and central normalized temperature L1 error is1.00153.
Removing this runaway does not establish general shock accuracy or convergence.

For a velocity element of width H containing a left dense sliver of fraction
delta, with negligible remaining density, its right-node diagonal mass is
M_rr≈rho_high*H*delta^3/3. The right pressure load scales as p_high*delta.
Thus a pressure load acts through very small represented inertia. Resolving
the material jump changes the spatial approximation rather than its force law,
mass, temperature, energy or arithmetic. This argument and the neutral result
do not alone prove the cause of the actual charged Ar9 trajectory.

The two preserved outputs are `material_density_step_buried`(26418 EXIT1) and
`material_density_step_resolved`(32480 EXIT0). Their frozen source is
`/private/tmp/frc-interface-resolution-20261017`, SHA
`597e00ce66a303deb1378a2c93b89f846454ccbed1ec7b7b2d0aa40cff086d2b`.
The independent whole-history readout and inspected figure are in
`outputs/material_density_step_interface_comparison`. Its companion documents
why the original prototype plots must not be used as matched complete endpoints.

## Actual reacting state and candidate

The captured final-stage Ar9 state has an8.39609-decade adjacent density jump
at±7.951239641mm, inside43.85654127um velocity elements. The dense material
occupies only3.068657% of each element. Other4.066812-decade jumps lie inside
740.413601um elements;7.470013-decade jumps already coincide with old nodes.

`density_velocity_refinement.refine_density_jumps` takes an explicit contrast
and adds corresponding existing fine faces to the mechanical mesh. It compares
log(N_relative/h); the physical number scale cancels. It never drops population,
old nodes or fine cuts and does not change physical mass/force/heat. The actual
acoustic bound still applies. There is no cutoff, spacing relaxation, hidden
retry or temperature floor. The contrast is a spatial-resolution parameter.

For the captured actual post-source initial state, contrast1e4 selects6 faces,
adds4 nodes and changes Ar9 from38 to42 velocity elements; minimum spacing
becomes1.345806842um. The old P1 field is interpolated exactly up to rounding.
The actual complete final5ns transport trial16298 completed successfully in
`outputs/reaction_density_jump_transport`, frozen root
`/private/tmp/frc-density-jump-20261017`, source SHA
`d06d0a56719348a33ec1a55eae352af3cc08ae9b3d1756197226672f12917eec`.
It retains all species/cuts and compares heat with T0*(h0/h1)^(gamma-1).
It omits the last source half by design and makes no coupled20ns claim.
The original full5ns control now completes241steps,2331.501s and Erel4.53608e-12.
Ar9 maximumT is370.160809eV, with369.748654eV maximum and29.730558eV volume-RMS
excess over reversible adiabatic compression. At its hottest cell, initialT
0.200055eV and width266.238nm become width90.0333nm and adiabaticT0.412156eV.
The enclosing velocity nodes accelerate from55.62/56.18 to54933.69/46477.30m/s.
Ar8 reaches14.176049eV before the last source return. This isolates the large
heat in transport, not the final reaction return; excess alone is not a shock
accuracy test. Complete states/steps/fields were fully decoded successfully.

Process56061 subsequently exits1 when disk storage fills while saving the next
level1 initial state, before that refined evolution begins. The original
summary and partial archive remain, with storage_failure.json distinguishing
the completed baseline from the incomplete comparison family. A fresh
continuation45403 loads the captured initial state without replaying the source.
It reveals that both midpoint rounds1/2 add ZERO nodes: actual Ar9 stays38,
and every initial array equals the completed control. The duplicate calculation
was retired with exit143 and explicit termination.json. No refined endpoint
or convergence claim. The density-interface trial is distinct: it adds4 nodes.

The completed control also has Ar5/6/7/8 peaks4.206/5.631/7.868/14.176eV.
New actual full5ns trial81727 therefore applies the same explicit1e4 density
criterion to all represented heavy species, including neutral Ar, while keeping
the32-element electron space. It loads the same captured initial state and
retains all N/U, old nodes and fine cuts. Actual Ar0..9 counts become
100/98/96/96/94/86/84/82/84/42. Maximum initial P1 representation discrepancy
is4.54747e-13m/s and measured relative energy change is0.0. This is still a
live stage trial, not an endpoint or coupled20ns claim. Output
`reaction_density_jump_heavy_transport`, frozen root
`/private/tmp/frc-density-jump-20261017`; same source SHA d06d0a56...17eec.
The unchanged original half-acoustic/plasma restriction and minmodQ apply.

A separate frozen full20ns argon trial applies the same explicit1e4 criterion
only to Ar9 after each source return, alongside existing heavy midpoint
enrichment. It uses the same variance projection,128-bit number/energy scaling,
physical masses,5ns source interval and all fine cuts. It completed successfully in
`outputs/reactive_material_argon_density_jumps_5ns`, frozen root
`/private/tmp/frc-source-density-jumps-20261017`, SHA
`57fc439319938d37ef9643d9b6a129f8dbf292061a121c9b8f4e87d493ff5b8b`.
Its complete20ns measurements and comparisons are at the top of this document.
The launch tool handle was lost in output truncation; terminal success is
established by the final summary/endpoints and process exit, without a separate
tool exit-code claim. Do not repeat this completed run without a concrete change.

The main optional integration is now evaluated; defaults remain unchanged.
Collect the all-heavy stage and continue source/space/threshold sensitivity
rather than treating one contrast as validated. Current accuracy,
molecular reactions/collisions, heat transport,
chamber boundaries and independent formation comparisons remain open.

## Independent mesh comparison complete

`material_density_step_cosine`(66102 EXIT0) clusters64/128/256 initial cells near the
internal jump and both vacuum boundaries while keeping the same physical
problem and duration. Uniform remains the default. All three runs finish the
full1.151043832us in3783/17733/79460steps. PeakcentralT rises
0.822741/0.909599/1.002660eV toward exact1.204242eV. Normalized temperature
L1 error decreases0.955088/0.784399/0.577695, density0.0123555/0.0062378/
0.0032400, velocity0.470952/0.352932/0.249352 and pressure0.0120684/0.0058522/
0.0030444. Finest maximumErel1.16029e-11. These L1 errors use initial scales,
not the exact field's L1 norm. Contact/shock location, shock temperature and
rarefaction tail remain underresolved; decreasing errors do not establish a
formal order or physical formation accuracy. Acoustic time steps change along
with spatial mesh. Original100000-step ceiling and duration were retained.
The final PNG was inspected. This run imported main and has now verified its
unchanged source closure; that live-source restriction has ended. Source SHA
`786801b144ccc5dcf9eb57dab4b3e1811e6ff385953f8d43bf57be040a6c4873`.

At05:14:56UTC the256-cell raw archive was still growing and only921714688bytes
remained free on the host. Its identified writerPID533/session66102 was
temporarily SIGSTOP-suspended to preserve memory and partial output before a
disk-full failure. `resource_pause.json` records this operational pause. No
physical state, duration or success flag changed. Resume with SIGCONT after
verifying identity and sufficient storage; do not restart or delete the partial
archive. Storage subsequently recovered and the exact process was resumed;
the complete family now succeeds. The host later filled during the original
Ar9 control's next initial export, producing the separate failure above. Other
argon processes were temporarily suspended and resumed with unchanged state.
A one-time clonefile operation additionally shares about620MB across92 identical
completed archives, preserving every path, byte and hash. Full cosine256 raw
history is2.1GiB; keep it local and avoid an additional Git copy while storage
is tight. Wall elapsed1928.343s includes export and a storage pause; it is not
a solver performance measurement. No scientific data or user poster was deleted.
