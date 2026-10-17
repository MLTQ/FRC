# Density interfaces and mechanical resolution

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
The actual complete final5ns transport trial16298 is running in
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
physical masses,5ns source interval and all fine cuts. It is running in
`outputs/reactive_material_argon_density_jumps_5ns`, frozen root
`/private/tmp/frc-source-density-jumps-20261017`, SHA
`57fc439319938d37ef9643d9b6a129f8dbf292061a121c9b8f4e87d493ff5b8b`.
At15ns it has completed the first three source intervals (367.368s elapsed).
At10ns it had completed325/203 transport steps and retained21161 density cells,
Erel1.71068e-11. Ar9 had30 velocity elements after the second return. This
prefix is not a completed physical result. Its launch tool handle was lost in
output truncation; processPID2415 was directly confirmed live. Do not restart.

No reaction-return default or main integration is changed by this candidate.
Collect the complete actual-stage and coupled results before adoption; then
measure threshold/resolution sensitivity rather than treating this one contrast
as validated. Current accuracy, molecular reactions/collisions, heat transport,
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
