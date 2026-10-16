# Spatial reacting-species state and electromagnetic fields

The spatial species state now has longitudinal charge/field feedback and a
transverse Maxwell extension. All physical species masses, three momenta and
total translational energies are retained. These components are not yet installed
at the moving formation solver's free boundary. Their fixed reflecting and
conducting intervals are verification geometries, not experimental chamber ends.

## Longitudinal Euler–Poisson result

Six complete nonlinear standing-wave histories are saved in
`outputs/atomic_electrostatic`, source
`8207447522e0740a74ecb14835b6042989489afcde0236cc920e6168998a16d2`.
The terminal summary reports success after 37.1613 s; every accepted state is
saved, and the scientific PNG was inspected. The launch handle was unavailable
after context recovery, so no process exit code is inferred from that summary.
No unit suite or energy/charge/population repair was used.

The independent adiabatic two-fluid dispersion relation and initial data are in
`examples/atomic_electrostatic.md`. Physical electron and argon-mass ion fluids
evolve on a 0.1 m interval at Nref=1e14/m³ and Ti/Te=0.2/2 eV, for 2.3 periods
(25.5235 ns). Both amplitudes 1e-3 and 1e-4 are exercised; the latter gives:

| Cells | Maximum E error / mode amplitude | Maximum J error / mode amplitude | Final mean density error / amplitude |
|---|---:|---:|---:|
| 16 | 0.1389895 | 0.1307185 | 0.1493319 |
| 32 | 0.0415455 | 0.0413549 | 0.0367140 |
| 64 | 0.00999065 | 0.0102031 | 0.00911662 |

Maximum timesteps are 70.90, 35.45 and 17.73 ps. Space and time refine together;
the continuum linear reference is not exact for finite-amplitude nonlinear Euler
evolution. Similar normalized errors at both amplitudes indicate mesh error
dominates here. The finest electric/current errors remain about 1%, not zero.

Maximum integrated nucleus/charge errors across all cases are 2.64e-16/3.33e-16
in the normalized integral units. The exactly absent neutral species stays zero.
Actual numerical continuity fluxes determine both the longitudinal Ampere current
and electric work. The sampled semidiscrete Ampere residual is at most 1.17e-15
of the mode rate, without imposing neutrality or correcting a field after its
update. These identities do not certify continuum accuracy by themselves.

For amplitude 1e-4, total particle plus face-lumped field-energy drift relative
to initial field energy is 1.07225e-4, 1.41065e-5 and 1.78656e-6. Using the exact
integral of the piecewise-linear Gauss field instead gives absolute budget drifts
2.89771e-9, 7.33853e-10 and 1.84067e-10 eV m per reference nucleus. This quadrature
difference decreases with refinement but is much larger than the matched-discrete
energy residual. Neither number is hidden or used to repair the state.

## Transverse Maxwell extension

`plasma_bench/atomic_maxwell.md` specifies the full spatial equations, staggered
field placement, source-work matching, RK4 update and physical light/plasma/gyro
timestep limits. `examples/atomic_maxwell.md` gives independent vacuum,
unmagnetized and magnetized circular-polarization references. Magnetic forces
and the induced longitudinal motion are included in the nonlinear evolution.
Fifteen complete histories finished in run86650 (exit 0, 92.079 s), preserved in
`outputs/atomic_maxwell`, source
`c2ab16c51682ca87da1ee3781a4e72d5717e52397cb1cb7e115277e7462450a3`.
Reflecting empty ghosts exposed an unnecessary 0/0 contact-state calculation in
the existing Euler flux, although the conservative fallback returned zero flux.
Both-empty faces now bypass that undefined intermediate calculation. No state
or nonempty face is changed. Full run99703 then exited 0 in 91.837 s, with all
15 histories in `outputs/atomic_maxwell_exact_empty`, source
`a8141a13a35d42704ec9a1e3d7766a09e0630bb819ab520771deb752a4826213`.
Every saved array agrees bitwise with the first complete run. Both scientific
PNGs were inspected; all fields, species states and PDFs are preserved.

Each trajectory covers 2.3 continuum periods. Maximum E/B vector errors relative
to the analytic mode amplitudes are:

| Mode | Cells | E error | B error | J error / mode amplitude |
|---|---:|---:|---:|---:|
| Vacuum | 16 | 0.0227685 | 0.0200719 | exact zero J |
| Vacuum | 32 | 0.00569152 | 0.00503985 | exact zero J |
| Vacuum | 64 | 0.00142292 | 0.00126122 | exact zero J |
| Unmagnetized | 16 | 0.0347937 | 0.0306890 | 0.0310040 |
| Unmagnetized | 32 | 0.00869782 | 0.00770209 | 0.00775313 |
| Unmagnetized | 64 | 0.00217434 | 0.00192734 | 0.00193835 |
| Circular, Bz=20 mT | 16 | 0.0331622 | 0.0350470 | 0.0321656 |
| Circular, Bz=20 mT | 32 | 0.00830295 | 0.00880310 | 0.00805710 |
| Circular, Bz=20 mT | 64 | 0.00207650 | 0.00220334 | 0.00201522 |

Plasma rows use E0=100 V/m. The additional E0=1000 V/m cases give essentially
the same normalized transverse errors, but their nonlinear longitudinal E is
100 times larger: the finest peaks are 0.0245121 V/m (unmagnetized) and
0.0368380 V/m (circular). That quadratic response is evolved, not suppressed.
Spatial/time refinement is simultaneous, so these are not independent temporal
and spatial order estimates. Their approximately fourfold error reduction is
consistent with second-order spatial accuracy.

The vacuum cases keep every species number, momentum and energy exactly zero,
while fields propagate at the physical light speed. This verifies empty-domain
Maxwell evolution, not a moving matter/vacuum interface. Across plasma cases,
maximum normalized integrated nucleus/charge drift is 7.63e-16/8.19e-16.
Finest matched-discrete energy drift is <=8.04e-10 of the initial transverse field
energy, but using piecewise-linear E integration gives about 4.02e-4 for the
vacuum and linear-polarization cases. The circular case's corresponding error is
4.66e-8. This difference is an explicit field-quadrature error, not proof of
continuum energy conservation to 1e-10. No budget is used as a state correction.

Time steps are at most 6.254, 3.127 and 1.564 ps at 16/32/64 cells. Resolving the
physical light/plasma scales this way is expensive and is not established as a
practical integrator for the much denser microsecond laboratory trajectories.
The completed homogeneous exponential-source run still shows substantial
reaction/field splitting error; see `ATOMIC_ELECTROMAGNETIC.md`.

## Remaining scope

The Eulerian matter flux can leak beyond a physical vacuum front; it must not
replace the tracked exact formation exterior. Spatial reactions/collisions,
Maxwell momentum/stress budgets, physical heat/photon transport, radial geometry,
open ends, real coil boundaries and a compatible moving interface remain to be
assembled and measured. A five-moment fluid wave is not kinetic validity,
formation, stability or confinement evidence. The original FRC goal stays open.
