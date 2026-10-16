# Full-vector moving plasma and current (2026-10-14)

The moving material solver now evolves all three axisymmetric components of B,
J and bulk velocity together, including material rotation, Hall/current inertia,
current stress, separate species entropy, prescribed resistance and heat. One
complete coarse trajectory finishes 7.121769 microseconds. This is a compact
plasma expansion with initialized closed flux, **not coil-driven FRC formation**.
The earlier finite-coil solver remains separate; penetrating vacuum fields have
not yet been coupled to this free material boundary.

## Equations and boundary contract

Use M=mi+me, mu=mi*me/M, w=vi-ve, J=enw, and delta=(mi-me)/M. The material
map carries poloidal displacement and the actual angle theta(R,Z). Electric
one-forms use G=F^-T G0; magnetic two-forms and curl G use F/det F. Thus angular
winding and changing cylindrical volume enter the same physical field/current
metrics as meridional motion. No angle or geometry derivative is frozen.

With Mn=integral dN G.T G, H=integral dV Bbasis.T Bbasis/mu0, and
D=-integral dV Bbasis.T curl(G)/mu0, the compatible equations retain

    e Mn w_coeff = -D.T b,
    H bdot = D E_material,
    Mn E_material = forcing + (mu/e) Mn wdot.

The forcing contains delta*w cross (B+mu/e curl w),
-(mu/e)*w cross curl v, species pressure work and eta*J. All moving Mn and D
derivatives remain in wdot. Physical vorticity includes 2*Omega+r*Omega_r;
the rotating coordinate frame does not remove this connection term.
Magnetic Maxwell stress and mu*n*w*w supply meridional force and angular torque.
Their work cancels the matching field/current geometry work. Bulk angular and
axial momenta are unforced. The full virial reference is I''=4E-2Wb.

Poloidal flux is represented by psi=r0²*a; reference Br=-psi_Z/r0 and
Bz=psi_R/r0. Both a and its normal derivative vanish at the outer side and both
ends. Together with the toroidal space this imposes **B=0 on the moving surface**;
exterior B, matter and conductivity are exactly zero. This restricted compact
boundary does not approximate a penetrating coil field or certify pointwise
carrier-current traces. No numerical atmosphere is present.

Species entropy uses the previously documented positive dual volumes and
first-order upwind numerical diffusion, with matching pressure work. Toroidal
scalar-gradient rows are exactly zero by axisymmetry. Constant eta=0.001 ohm-m
is prescribed; actual eta*J² heats the species with mass-weighted partition.
Conditional Maxwellian exchange uses their evolving temperatures. There are no
entropy/density/log floors, omitted edge heating or energy corrections.

## Complete trajectory and independent comparison

The run uses R=35 mm, half-length=70 mm, central density1e18/m³, physical Ar+
and electron masses, n=nc*(hr*hz)^1.5 and initial Te=Ti=hr*hz eV, where
hr=1-r²/R² and hz=1-z²/L². Initially there is no bulk motion. The analytic fields
come from psi=(Bp/2)*r²*hr⁴*hz⁴ and Btheta=Bt*r/R*hr³*hz³, with Bp=200 microtesla
and Bt=50 microtesla. Analytically differentiated curl(B)/mu0 is supplied to the
physical magnetic/carrier projection. The independent initial force is
-grad(p)+J cross B-div(mu*n*w*w), with all cylindrical stress connections.

| Measured quantity | 4x8, degree3, q6, 2000 steps |
|---|---:|
| Initial B / J relative volume L2 errors | 19.617% / 35.856% |
| Initial acceleration RMS error / (C/R), C=10e/M | 0.0100964 |
| Maximum relative total-energy change | 2.32664e-6 |
| Maximum independent virial error / initial inertia | 2.00949e-7 |
| Maximum absolute total axial / angular momentum | 2.61e-24 kg m/s / 7.40e-30 kg m²/s |
| Largest generated toroidal speed / final | 8.29662 / 7.31690 m/s |
| Observed Joule heat | 434.214 nJ |
| Final peak B / J | 0.723917 microtesla / 81.4865 A/m² |
| Final particle-mean Te / Ti | 0.375864 / 0.367898 eV |
| Minimum species entropy over entire run | initial value, -27.6310211 |
| Largest outflow rate times dt | 0.0314989 |

The initialized poloidal field decays strongly; final peak Br/Bz are
1.35e-11/4.21e-11 T while peak Btheta is7.24e-7 T. All are evolved values.
This behavior is not a demonstrated converged decay rate. Energy conservation
does not erase the large initial current error, upwind diffusion or uncertain
cold-edge transport. No laboratory design, confinement or 3-D stability follows.

## Evidence and next simulation

`outputs/material_vector_free_r4_t2000` preserves every accepted state and frozen
source closure5744fb58115ea903ba9a683959e4e47b30b331d75981d11bb672ae7e4fa743e6.
Scientific handle20443 exited0; wall time441.26 seconds, peak memory142.2 MB.
The [interactive display](../outputs/material_vector_display/index.html) contains
five actual states, linked 3-D magnetic paths/current samples, signed components,
shared/per-state scales and actual moving boundary. PNG/PDF are saved alongside.
Paths use spatial interpolation, never particle tracking or topology acceptance.
Browser time/component/scale/playback/rotation/zoom/reset controls were exercised;
no warning/error logs were returned. The scientific PNG was visually inspected.

A full 6x12,degree5,q10,2000-step run is LIVE as handle89864 in
`outputs/material_vector_free_r6_p5_q10_t2000`. Its independent **initial** B/J
errors fall to0.4320%/0.8437%, and acceleration error to0.0018903 of C/R.
Completion and full-trajectory accuracy are pending. Poll the same handle;
do not restart it on an observation timeout or edit its live source closure.
After the current run, use its actual physical comparison to direct the next
solver change, then couple penetrating moving vacuum/finite coils. No new
verification framework or unit tests accompany this extension.

## Completed higher-order trajectory (collected2026-10-15)

The6x12,degree5,q10 trajectory completed all2000steps/7.121769us in11483.19s.
The old process handle89864 expired overnight, but its saved terminal summary
confirms success. No restart was made. All raw states and the original frozen
source5744fb58115ea903ba9a683959e4e47b30b331d75981d11bb672ae7e4fa743e6
remain in outputs/material_vector_free_r6_p5_q10_t2000.

Initial B/J errors0.4320%/0.8437% and force RMS error0.00189028 of C/R compare
against the independent analytic initial fields/force above. Maximum full-run
energy change1.49326e-5 and independent virial error1.74261e-6 remain measured,
not repaired. Maximum absolute angular/axial momentum7.75e-30kg m²/s and
1.36e-24kg m/s. Final Joule heat479.663nJ, peak B0.675541uT, peak J102.315A/m²,
meanTe/Ti0.376493/0.367776eV. Final peak toroidal velocity18.1868m/s, versus
7.31690m/s on the older coarse run, and Joule heat is10.47% larger. Those changes
preclude claiming converged rotation/current from integrated conservation.
Both resolution and field degree changed; this is not a measured convergence
order. Compact initialized flux remains distinct from coil-driven formation.
