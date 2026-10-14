# Configured metre-scale biased-coil experiment (2026-10-14)

The simulation now accepts an explicit JSON experiment: tube/plasma geometry,
gas inventory, initial species temperatures, bias, serialized pulse shape, mesh,
prescribed resistance and duration. The supplied case is
`configs/material_tube_bias_pulse.json`. Run it with
`python -m examples.material_experiment CONFIG OUTPUT`.

## Initial assumptions and actual evolving model

The example uses a1m-long,15cm-diameter nonconducting/nonmagnetic tube, matching
one of the user's provisional geometries. Gas is Ar at0.1mTorr/300K before an
ASSUMED preparation step. Its5.68815e16 atoms are all assumed present as a compact
singly ionized plasma, R35mm and half-length350mm, with n=nc*f^(3/2) and
Te/Ti=(10,1)*f eV, f=(1-r²/R²)(1-z²/L²). Actual volume normalization gives
nc=8.96143e19/m³. The neutral fill density3.21883e18/m³ is not that peak density.
Preparation, dissociation/ionization, neutral losses and their energy are absent;
this is an exploratory initial state, not a claim that this preparation works.

Two prescribed circular loops at r100mm,z±50mm supply a100uT centre bias. A
-1mT pulse starts at0.25us, rises1us, holds0.5us and falls1us. The full3.5us run
includes return to the bias. All loops share the one current. Bias is highly
nonuniform along the long tube: initial plasma quadrature Bz ranges3.157uT to
102.397uT. Helmholtz spacing does not make a metre-long uniform-field region.
The amplitudes are abstract simulation fields, not hardware recommendations.

The actual full material/current/species/magnetic-vacuum equations evolve with
physical masses, Hall/current inertia and stress, free side/end motion, rotation,
Joule heat and conditional Coulomb temperature exchange. eta=1e-4ohm-m is an
explicit prescribed scalar. All exterior matter/conductivity stays exactly zero.
The simulation stops on tube contact at any stage; it does not invent wall/sheath
forces. The inherited upwind entropy diffusion and local current time sensitivity
remain. Initial minimum temperatures are0.001416eV/0.0001416eV, minimum exchange
Coulomb log3.744: even the initial cold edge is not a validated physical plasma.

## Current spatial accuracy

The completed6x12 coil pulse6211 refines the preceding4x8 fixture with the same
physical inputs and2000steps. Final Joule heat changes5.49395 ->5.44027nJ. Actual
volume-weighted field/current samples at common material labels give:

| Time, us | Total B difference | Total J difference | Velocity difference |
|---|---:|---:|---:|
|1.78044|0.262%|30.26%|0.724%|
|3.56088|0.628%|8.56%|0.747%|
|7.12177|23.62%|60.38%|0.796%|

Differences use fine-run actual volumes and common10-point Gaussian quadrature on
the union of reference cells. Maximum final physical-position mismatch is0.319mm.
Final azimuthal-current difference59.43% remains substantial; dominant source B
can hide induced-field/current error during the pulse. These two grids do not
establish convergence. Independent peak-pulse exterior coil linkage discrepancy
improves0.3791% ->0.3074%. The inspected plot/data is in
`outputs/material_coil_spatial_comparison_q10`; actual field basis reconstruction
uses `material_plasma_samples`, without interpolating a quadrature cloud.

The8x16 version of this same full pulse is running to determine the trend in the
actual current distribution. This is a concrete simulation accuracy problem;
no new test suite, recursive verification framework or pass/fail threshold.

## Experimental baseline to pursue

Peng et al., *Formation studies of field-reversed configurations on the HFRC-F
device*, report a typical density4e20/m³ and75us lifetime with bias-15mT,
pre-ionization15mT and main35mT. Their internal Bz and line-integrated density
histories are useful experimental targets. The paper reports that most ionization
occurs after main-field triggering and that preparation timing changes trapped
flux. It therefore cannot validate the present fully pre-ionized3.5us case just
by matching a field amplitude. Actual geometry, input history and preparation
must be matched first; no digitized measurement or model-to-experiment fit is
claimed. [Primary conference paper](https://indico.fusenet.eu/event/28/contributions/101/attachments/115/999/Formation%20studies%20of%20field-reversed%20configurations%20on%20the%20HFRC-F%20device.pdf).

Sgro et al. find that measured formation flux loss and temperatures require
transport and ion-heating treatment beyond a classical resistivity assumption
in their experiment. That is evidence to investigate transport and heat partition,
not permission to tune coefficients until this model agrees. The accessible
abstract also identifies significant current in tenuous plasma outside the FRC;
a real modeled plasma region must be distinguished from an invented numerical
atmosphere. [Primary article abstract](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.26.3564).

## Complete first trajectory and current continuation

The first metre-scale trajectory20691 completed2000steps/3.5us successfully in
1141.29s, peakRSS212.39MB. Frozen source identifier:
`e171243585fdd4bf3f8a856650aa4b279d2a0d72505073a556d7eb4244b205b7`.
Output: `outputs/material_tube_bias_pulse_r4_z16_t2000`.

Peak current components(r,theta,z) are(271.94,84732.49,2609.58)A/m². Final
plasma radius54.350mm remains below the75mm tube wall; half-length352.340mm
remains below500mm. Final mean Te/Ti=3.49554/0.375728eV, Joule heat195.603uJ,
source work156.495uJ, induced-field energy22.1122uJ and current energy0.568513uJ.
Initial plasma energy is89.5071mJ. Maximum external and mutual/source work
residuals relative to initial energy are3.8827e-7 and4.1213e-7. Their smallness
is no local-current, transport or formation certificate.

At3.5us the imposed centre field has returned to+100uT while actual plasma Bz
ranges-142.559 to+171.436uT.77.22% of sampled plasma volume has Bz opposing the
local source Bz. The displayed flux contours and field lines show the coarse
solution's structure; no resolved separatrix, stability or confinement is claimed.
The maximum observed electron-drift/thermal-speed ratio reaches5.9535; the final
cold-edge sample is still1.4946, with ion-ion Coulomb log0.01593. These diagnostic
results contradict treating a Maxwellian collisional closure there as established.
No neutral/ionization or kinetic edge model has been supplied.

The completed actual-state3-D view and signed Bz/Jtheta/flux sections are in
`outputs/material_experiment_display`. They include the actual exterior field,
coil locations and tube outline. Current remains confined to plasma.12 actual
saved frames cover bias, pulse transitions and relaxation; color scaling is
explicit. 3-D lengths share one physical scale;2-D sections explicitly expand
the radial scale for readability. Display38726 exited0; PNG/PDF figures inspected. BrowserTab8 at127.0.0.1:8790
(server89915) is marked deliverable. Playback, all current components, initial
zero current, scaling, rotation/zoom/reset were exercised without browser errors.
`browser_view.jpg` preserves the inspected final view.

## Pulse/control result and stronger-field continuation

Bias-only control63566 exited0 after all2000steps/3.5us (1150.89s). Its final
radius54.261733mm differs from the1mT pulse's54.350118mm by only88.39um.
The control has no sampled axial-field reversal: Bz2.751..123.465uT, versus
-142.559..171.436uT in the pulse. Final Joule heat is0.520183uJ versus195.603493uJ;
peak final Jtheta3.53771 versus16.51733kA/m². The pulse produces magnetic response
and heating but does not arrest expansion. Both mean ion temperatures remain
about0.376eV. The actual signed field/flux, current, heat and radius comparison is
`outputs/material_tube_pulse_control_comparison/pulse_vs_bias.png` and PDF.
The complete successful trajectories, source hashes and figure were inspected.
This is a physical control, not a new verification framework.

The third weak-pulse grid42299 exited0:8x16,p3,q6,all2000steps,2283.64s.
Actual6x12-to8x16 current differences are16.73%,8.94%,48.63% at quarter/half/end;
final azimuthal current differs47.28%, total B10.05%, velocity0.3684%.
The4x8-to6x12 final current difference was60.38%. This is an improving final
trend but not local convergence; intermediate-time differences are not uniformly
monotone. Peak-pulse exterior/volume coil linkage differs0.29949%, versus0.30738%
on6x12. Existing exact-basis comparison: `outputs/material_coil_spatial_r6_r8_q10`.
No solver states were fitted to either energy or these observations.

An independent initial pressure scale helps choose the next physical pulse:
p_center=nc*e*(Te+Ti)=157.936Pa and sqrt(2*mu0*p)=19.923mT. This is only a static
magnetic-pressure comparison; it does not predict dynamic capture or required
hardware. The1mT pulse is weak compared with this initially assumed pressure.
`configs/material_tube_bias_pulse_20mt.json` keeps every initial plasma input
and waveform time unchanged, raising the opposing prescribed pulse to20mT.
The20mT run52842 exited0 after all2000steps/3.5us in1197.04s, peakRSS202.95MB:
`outputs/material_tube_bias_pulse_20mt_r4_z16_t2000`, same frozen source e1712435.
The minimum recorded midplane radius is31.3557mm at2.10us (initial35mm). It then
rebounds to53.1791mm; maximum boundary radius reaches68.6460mm, compared with
54.2617mm in the bias-only control. Thus central compression is transient and
outer expansion is greater. Half-length ends352.3805mm. Final Joule heat85.5213mJ
and external work75.3002mJ compare with89.5071mJ initial plasma energy. Mean
Te/Ti ends7.12060/.439237eV. Final plasma Bz spans-5.28416 to+.929770mT while
source bias is back to+.100mT. No sustained confinement or experimental FRC claim.
The surface maximum is not a bulk confinement measure. Actual particle-weighted
r_rms=sqrt(sum(dN*r²)/sum(dN)), using Q1-mapped quadrature positions and fixed
particle counts, starts18.7082mm. At3.5us it is27.4959mm (bias only),27.4887mm
(1mT),25.6110mm (20mT). The stronger pulse therefore reduces final bulk radial
extent6.85% relative to control while still expanding36.90% from initial.
Axial RMS extent goes142.8868 ->143.7158mm (20mT), versus143.4984mm control.
These actual moment integrals are in `material_tube_20mt_control_comparison/particle_extent.json`.
No escaping-particle or confinement-time model is present.

Maximum external/source-work residuals/E0 are2.72010e-5/4.95757e-5. Peak electron
drift/thermal-speed ratio reaches116.93 and ends28.58; the Maxwellian collisional
edge closure is not validly established. The larger pulse has not been time or
spatially converged. Its displayed morphology is an exploratory solver result.
The actual comparison figure is `outputs/material_tube_20mt_control_comparison`.
Dashed radius follows the material midplane; solid radius is the largest boundary
value, so local compression cannot conceal outer expansion. Both compared field
sections share a signed scale. No modified heat, atmosphere or state correction.

The completed12-frame3-D view is `outputs/material_experiment_20mt_display` and
http://127.0.0.1:8790/20mt/ (existing server89915, relative symlink in the older
1mT view; original data remains unchanged). BrowserTab8 now shows this run and is
marked deliverable. New live radius labels distinguish middle/max geometry.
Time selection and signed azimuthal-current selection were exercised, then the
final magnitude view restored. No browser warnings/errors. PNG/PDF scientific
figures and browser_view.jpg were inspected; the screenshot preserves the view.

A candidate time integration in `material_split_current` separates the existing
upwind numerical diffusion and treats its discontinuous matched pressure force
at zero current implicitly. Its exact equations and numerical diffusion are in
the companion. Actual Hall, inertia, stress, species, heat and vacuum terms remain.
The first two full-pulse attempts failed early (both at step2);
all accepted states and frozen sources remain in `material_coil_split_r4_t2000`
and `_b`. The corrected active free-row solve is evaluated only through complete
pulses: The c runs ended at steps1656/808 after a newly hit bound repeatedly released
caused a zero-length solve step. Captured actual matrix: `material_split_failure_diagnosis`.
The retained working-set solve reduces its equation residual to4.34e-19 without
changing equations or bounds. Both corrected runs14252/45710 have now completed (full results below). No adoption or temporal
accuracy is claimed until the full saved currents can be compared. Direct
-Idot*lambda work remains independent; ambiguous instantaneous source-response
work is absent for this nonsmooth method rather than fabricated.

LIVE89864 remains `material_vector_free_r6_p5_q10_t2000`, the older full-vector
high-order trajectory. Preserve every live source closure and poll exact handles.
No unit tests. The full FRC goal, physical transport/ionization and confinement
remain open. Next use complete stronger-pulse motion/field results to choose
physical simulation improvements, with the current comparison informing accuracy.

## First completed implicit upstream pulse

Corrected run45710 completed all1000steps/7.121769us (539.42s,176.24MB RSS),
`outputs/material_coil_split_r4_t1000_d`, source
`d56d01a8aa762d4b58f241fb8abd3b5636278de928764c24497fdca1de1a0eeb`.
Final J(r,theta,z) peaks are.704557/9.91443/1.36981A/m²; heat5.49399994nJ,
external work5.53647531nJ. Maximum external-work residual/E0 is1.06014e-9;
maximum observed diffusion energy-exchange residual is1.58602e-20J. The species
minimum entropy falls below its initial value only by2.13e-14/1.60e-13, consistent
with roundoff here; no clipping or lower-bound reset occurred.

Against the completed original RK4 run at the same1000steps, final total B/J
change.9814%/1.8543%; Jr/Jz differ19.80%/18.11%, Jtheta.1504%. Quarter/half total J
differences are.09417%/.03050%. This shows that weak late currents depend on the
time method; it does not prove which current is accurate. Inspected plot/data:
`outputs/material_coil_split_vs_rk4_t1000`. Source-response integral stays absent
for the nonsmooth method; direct waveform work is compared.

The2000-step version14252 has now completed with the same corrected source;
its full time-halving comparison follows below. LIVE89864 remains the older high-order vector trajectory.
The new step stays opt-in pending that result. No additional verification suite.

## Completed time-halved upstream pulse

Run14252 also exits0: all2000steps/7.121769us,1066.17s,172.00MB RSS, same d56d01a8
source. The complete1000/2000-step physical comparison is
`outputs/material_coil_split_time_comparison`; scientific PNG/PDF inspected.
At quarter/half/final duration, total B differences are.0002917%/.0005029%/.21196%,
total J differences.015047%/.009147%/.38724%. Final Jr/Jtheta/Jz differences are
3.7249%/.018683%/4.0554%, compared with15.28%/.150%/14.99% in the original RK4
pair. Final Btheta difference remains6.0173%. The new integration substantially
reduces the measured late-current time sensitivity, but does not erase it.
Maximum position difference is1.246nm and final velocity difference3.03e-8.

Fine-run maximum direct-work residual/E0 is2.65160e-10, about fourfold smaller
than the coarse1.06014e-9. The independently trapezoidal work observation also
contributes to this residual, so this ratio alone is not a temporal-order proof.
Fine final heat5.49397388nJ, peak J(r,theta,z)=.689890/9.913262/1.335545A/m².
Maximum observed diffusion exchange residual1.72340e-20J; species minimum entropy
changes relative to initial are no lower than-2.13e-14/-2.91e-13. No clipping.

Keep the new method opt-in: this evidence is the weak10uT pulse, not the stronger
configured experiment or its nonconvex pressure branches. Spatial current errors,
physical edge/transport and ionization remain open. All new scientific jobs are
terminal; only the older high-order full-vector89864 remains live. Next prioritize
physical simulations with independently prescribed bias/formation coil groups
and defensible transport/preparation, evaluating concrete changes end to end.
Do not turn the remaining error into an exhaustive replay/reporting milestone.
