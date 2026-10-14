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

LIVE63566: `outputs/material_tube_bias_only_r4_z16_t2000`, same assumed inventory,
geometry and initial state with zero pulse amplitude, for a complete bias-only
physical control. LIVE42299: `outputs/material_coil_free_r8_t2000`, the next
spatial refinement of the earlier complete weak pulse. LIVE89864 remains the
older high-order compact full-vector trajectory. Poll the exact handles and
preserve their source closures; never restart because an observation timed out.
Next compare actual pulse/control field, current, expansion and heat; collect
spatial trend, then address numerical/local-current and physical transport limits.
The full FRC goal remains open. No unit tests or new verification framework.
