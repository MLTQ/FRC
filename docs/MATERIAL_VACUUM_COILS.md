# Moving plasma, penetrating vacuum and prescribed finite coils (2026-10-14)

The free material solver now has a magnetic-only exterior following its actual
side/end deformation. A complete unforced plasma/current/species trajectory
finishes with penetrating poloidal flux. Finite circular-coil forcing is also
implemented, with two complete-pulse trajectories finished. These are
coarse component trajectories, not converged FRC formation or laboratory inputs.
The user's metre-scale tube, selectable gas and neutral fill remain exploration
requirements; the centimetre-scale fully ionized fixture below does not validate
them. No unit tests or new verification framework were added.

## Actual exterior and current coupling

Five magnetic-only patches extend the moving material side and both ends to a
fixed rectangle, with continuous analytic multipoles and energy to infinity.
Q1 boundary displacement extends to zero at that outer box. This coordinate rule
adds no exterior mass, pressure, conductivity or current. Folded geometry rejects.
Flux, its gradient and the SAME energy use the actual moving map at every stage.

Poloidal coefficients are nodal psi in x=r0²/2,z0 with psi(axis)=0. The exterior
uses the identical polynomial trace everywhere, including the ends. This is a
distinct space from the earlier compact r0²*a interpolation; coincident nodal
values alone would not match their traces. Btheta has free reference boundary
coefficients, allowing them to compensate winding of penetrating poloidal flux.
Exterior Btheta is zero; tangential matching is a weak condition whose finite
errors remain to be measured, not an assumed resolved surface-current layer.

The full magnetic metric is H=Hinside+T.T*Hvac*T. Vacuum internal coordinates
minimize magnetic energy at actual boundary flux; Hvac_dot differentiates that
same geometry, not a frozen exterior. For toroidal electric coefficients,
D[:,Etheta]=-H[:,poloidal] and Ddot=-Hdot; other columns use the interior curl
work. Thus weak Ampere, Faraday and finite current inertia include exterior flux
work. Adding exterior energy to diagnostics alone would omit this coupling.

Minus the exact minimized vacuum-energy geometry derivative supplies the exterior
force on actual material coordinates. It cancels the corresponding metric work,
without a residual-fitted force or state repair. Full physical masses,
Hall/self-drift, current stress, rotation, species entropy/pressure and heat remain.
The existing upwind numerical diffusion and its measured time sensitivity are
inherited. Neither exact vacuum nor discrete power pairing establishes cold-edge
transport validity or pointwise current traces.

## Independent moving-vacuum trajectory

A fixed physical analytic dipole is sampled outside a prescribed nonaffinely
expanding cylinder: midplane radius grows25%, poles extend15%, corners stay fixed.
Separate physical boundary integrals give exact dipole exterior energy, Maxwell
traction and moving-domain energy rate. Three complete21-state trajectories
compare actual solved fields/forces against those references:

| Maximum error over trajectory | 4x8 / vacuum4 | 8x16 / vacuum8 | 16x32 / vacuum16 |
|---|---:|---:|---:|
| Exterior field, relative volume L2 | 9.059% | 1.379% | 0.2659% |
| Nodal Maxwell traction, relative L2 | 4.599% | 0.5417% | 0.08055% |
| Exterior energy, relative | 1.778% | 0.2037% | 0.02940% |
| Moving-domain power, relative | 3.553% | 0.3466% | 0.01341% |

All use degree3/q6; multipole counts also increase6/10/14, so this is joint
refinement, not an isolated polynomial-order claim. The independent five-point
physical sampler (including axis, intended coil locations and beyond the outer
box) has final field errors4.485%,0.6093%,0.1621%; flux errors0.458%,0.04665%,
0.002783%. It inverts the actual element geometry and evaluates its field basis,
without interpolating a quadrature cloud. Finite derivative jumps at element
interfaces are not averaged away. The inspected PNG/PDF and values are in
`outputs/material_vacuum_motion_comparison`.

## Complete freely moving plasma

`outputs/material_vacuum_free_r4_t2000_corrected` completes2000 steps/7.121769us.
It uses R35mm,Lhalf70mm,nc1e18/m³,physical Ar+/electron masses and the existing
compact density/temperature profile. An analytic spherical volume current within
a=.7R supplies an initial200uT central field and a dipole exterior. Its psi and
first two derivatives match at the current-support surface; J vanishes there,
so the continuum initial field has no singular current sheet. This is initialized
flux, not a formation experiment.

Initial B/J errors are19.75%/14.31%; initial acceleration RMS error is0.005335
of10e/(M R), and vacuum-energy error is-0.1544%. Over the actual free trajectory,
maximum total-energy change is3.0043e-6 and all-space virial error3.4256e-7.
Maximum paired power residual is1.63e-13 W; these balances do not erase spatial
errors. Final Joule heat336.62nJ, mean Te/Ti0.37379/0.36836eV, peak velocities
(r,theta,z)=(3732.47,0.5934,2213.05)m/s. Vacuum geometry stays positive, with
minimum final determinant0.81346. Cold-edge/transport validity remains unresolved.

Scientific handle40493 exited0. Its frozen source closure is
05d7546a0e3b30d0b0bc586b04157f16ff86a9615b169b82ad381f11f80c936e.
The run archive used the temporary entry-point name material_vacuum_expansion;
the new working entry point is material_penetrating_expansion. The original
nonmagnetic affine driver at the old name has been restored unchanged. Two
initialization failures (metadata key and a name collision with the inherited
vacuum response matrix) remain preserved; neither produced a physical trajectory.

## Complete finite-coil pulses

Analytic circular loops remain separate prescribed fields: total B=Bind+I*Bunit,
source material E=-Idot*Aunit+I*v cross Bunit. Both this electric drive and the
actual J cross source-B force/torque enter the complete coupled equations.
Every RK stage uses its actual source time. Filaments crossing/touching the
moving material reject. The physical mixed-row tolerance remains1e-10.

Reciprocal linkage lambda=integral J.Aunit uses actual current. Its derivative
includes current rate, material electric-map rate and the exact spatial
derivatives of Aunit. The two physical response budgets are

    d(Efluid+Kcurrent+Wbind)/dt = -Idot*lambda,
    d(Efluid+Kcurrent+Wbind+I*lambda)/dt = I*lambdadot.

Only the divergent prescribed filament self-energy baseline is excluded. Work
is computed from source/linkage, never fitted to the observed energy change.
Exterior induced flux at physical coil locations gives an independent linkage
comparison; finite-space disagreement must remain visible, not repair the work.

Both 1000/2000-step runs (30464/26556) completed successfully. They use r70mm,
z±35mm loops, a10uT nominal vacuum-center pulse (1.5us rise,1us hold,2us fall),
followed through7.121769us. This is no DC-bias/formation study. Prescribed eta is
0.001ohm-m, not inferred physical transport. Actual particles are the same
centimetre-scale fully ionized fixture above, initially unmagnetized and at rest.

| Measured quantity | 1000 steps | 2000 steps |
|---|---:|---:|
| Maximum external-work residual / initial energy | 7.5990e-9 | 3.5432e-9 |
| Maximum mutual/source-work residual / initial energy | 1.1236e-8 | 3.9737e-9 |
| Peak azimuthal current, A/m² | 136.243 | 136.095 |
| Final Joule heat, nJ | 5.4939498 | 5.4939474 |
| Final source work, nJ | 5.5364474 | 5.5364290 |

Both preserve the initial global entropy minimum. At1.780us, same-material-label
relative L2 time differences are0.00193% in total B and0.1013% in total J. At the
end, those differences are0.624% and1.498%; however, individual current components
(r,theta,z) differ15.28%,0.150%,14.99%, and Btheta differs16.24%. Dominant-component
agreement hides weaker late-time dynamics. Final bulk-velocity difference is
7.87e-8 relative and maximum nodal-position difference2.31nm. These results do
not establish complete time convergence, much less spatial or physical validity.

At2.001us, the volume-current linkage is-1.9266152nWb while independent physical
vacuum flux at both coils gives-1.9339181nWb, a0.3791% discrepancy. This is visible
finite-space reciprocity error, not a correction to either source-work budget.
The inspected actual waveform/current/heat/work and Bz/Jtheta maps are in
`outputs/material_coil_comparison/material_coil_pulse.png` and its PDF;
`comparison.json` records source identifiers and componentwise time differences.

LIVE6211: `outputs/material_coil_free_r6_t2000` repeats the complete pulse on6x12,
p3/q6, with the same vacuum4/modes6 and2000steps. Poll its exact handle and compare
physical current/field/motion and coil linkage after completion. The older
compact6x12,p5,q10 full-vector run89864 also remains live and untouched. Do not
restart from observation timeouts or edit either source closure.
The completed upwind time control56380 is discussed in MATERIAL_THERMAL_COUPLING.md:
its final current remains14.36% time sensitive. Upstream entropy-pressure switching
is a candidate issue to investigate, not a reason to discard those measured errors.
Keep the overarching goal and frc-vacuum.3 open; local current/time accuracy,
physical transport/ionization, formation and experimental controls remain.
