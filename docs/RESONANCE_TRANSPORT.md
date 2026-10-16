# Resonance photons in the provisional neutral-argon tube

## Complete transport histories (2026-10-16)

Run 92249 exited 0 after 524.793 s. `outputs/resonance_transport` contains 20
complete cases, per-packet exits, first-eight-packet event histories, population
histories, the source/input archive and inspected PNG/PDF. Source SHA256:
`c896aed3898d0b2e76bab5f09f235ac64326d9ad173a2f72d2a84d9f730d84bc`.
NIST radiative input: `bf04549cd6dae143ff65ee92700902f1eea0567003b005fea5310dd598d5f2db`.

The physical component follows a weak, impulsive centre excitation through a
stationary, uniform 300 K neutral cylinder, radius 0.075 m and length 1 m.
All photons eventually cross a side or end; exterior opacity is exactly zero.
This is a radiation Green function, not a prepared discharge or formation run.
No timestep, spatial grid, escape-factor fit, packet killing or normalization is
used. Output times only sample exact emission/flight/absorption intervals.

Absorption is Einstein line strength times a thermal/natural Voigt profile.
Default redistribution samples the absorbing atom's conditional thermal
velocity and preserves photon frequency in that atom's rest frame. An explicit
complete-randomization comparison independently redraws each emission frequency.
Read the module companions for the conditional density, boundaries and exclusions.

## Measured escape

4096 packets per case. Times below are sample means in microseconds; uncertainty
is one sample standard error. Percentages describe radiation still inside the
gas at 3.5 us, including quanta stored as atomic excitation.

| Pressure equivalent | 104.8 nm correlated mean | 106.7 nm correlated mean | Retained at 3.5 us, 104.8 / 106.7 nm |
|---|---:|---:|---:|
| 0.1 mTorr | 0.12614 ± 0.00168 | 0.11808 ± 0.00157 | 0 / 0% |
| 1 mTorr | 1.40649 ± 0.01892 | 1.61268 ± 0.02116 | 6.64 / 8.86% |
| 10 mTorr | 9.43456 ± 0.12912 | 18.21295 ± 0.24654 | 73.27 / 86.77% |
| 100 mTorr | 70.46226 ± 0.95794 | 116.01103 ± 1.61133 | 96.46 / 97.61% |

At 10 mTorr complete randomization instead gives means 1.26921 / 7.24598 us and
retained fractions 5.74 / 64.06%. The frequency correlation assumption matters.
This comparison does not represent a collision-broadened gas: both modes retain
the same natural/thermal widths and exclude collisional destruction.

At 106.7 nm, 10 mTorr, independent seeds with 1024 / 4096 / 16384 packets give
means 18.32431 ± 0.47814 / 18.21295 ± 0.24654 / 17.68093 ± 0.11805 us.
The last two differ by 1.946 combined standard errors; their retained fractions
are 86.62 / 86.77 / 86.61%. These are statistical sensitivity measurements,
not deterministic timestep convergence. Do not select the closest seed.

## Independent references and budgets

Vacuum means are 2.24812 ± 0.02995 ns and 7.80342 ± 0.11912 ns, versus exact
natural-wait-plus-ray-flight means 2.25398 and 7.95004 ns. Maximum differences
from the analytic escape CDF are 0.00948 / 0.01430. For all 18 absorbing cases,
first-flight escape differs from independent frequency/angle Beer-Lambert
quadrature by at most 1.913 binomial standard errors. This checks free flights
and boundaries, not the physical accuracy of many-scattering redistribution.

Every run finishes every packet (maximum 332558 absorptions; latest exit
753.310 us). Excited + flying + escaped packet imbalance is exactly zero at all
sampled times, with no negative occupancy. Each carries constant h nu0: this is
an inventory budget, not full plasma/atom/photon energy conservation. Recoil,
frequency-shift energy exchange, atom motion, polarization, collisional mixing
and quenching, branching and surface reflection/transmission remain absent.
An exponential natural holding time even in far wings is a declared semiclassical
approximation. None of this source physics is coupled into the moving plasma yet.

The unbounded Voigt wings expose a further limitation of the complete-randomization
comparison. For 104.8 nm at 100 mTorr, 73/4096 exiting packets exceed relative
detuning 0.001, 7 exceed 0.01 and one reaches 0.22337. At 10 mTorr, 4/4096
exceed 0.01 (maximum 0.09345). These rare draws are retained, not clipped, and
cannot be called physically faithful narrow-line photons. No correlated case
exceeds 0.001 in these finite samples; that is not a bound on the formal tails.

## Experimental comparison to build next

ORNL-TM-4368, *Transport of Resonance Radiation in Argon* (1973), supplies
1048 A time-dependent measurements. The original report and visually checked
Table 1 are saved under `data/atomic/radiation_reference`. Its 7 cm × 2.22 cm
cell, axial excitation and narrow side view differ from this centre-pulse tube.
The published early-time detector decay rate is not inverse global mean escape
time. Reproduce the source and detector, preserve fit-window dependence, and
retain the paper's impurity flags before claiming experimental agreement.
Quenching, pressure broadening and detector spectral acceptance also need review.
Keep the atomic/moving coupling and full FRC goal open. No unit tests were run.
