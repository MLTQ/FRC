# Argon radiation: finite-cell detector comparison

## Completed calculations (2026-10-16)

The radiation component now supports distributed excitation, a narrow side
view, and sourced resonant collision broadening/decorrelation. This is a weak
impulsive source in a stationary neutral bath, not an evolving FRC plasma.
Geometry matches known ORNL-TM-4368 dimensions: length 7 cm, radius 1.11 cm,
side aperture diameter 0.714 cm and axial slit width 50 um. The assumed uniform
3 mm beam disc is based on the upstream collimator, not a measured beam profile.

Direct-position study 27479 exited 0 after 294.751 s: 16 absorbing batches and
a 262144-history vacuum trajectory, source
`67932429ed2f3e22ebbdc8e702bb861fe9e739162301288d3f1978b15e6fccee`.
`outputs/resonance_experiment` preserves all histories. It uses 300 K with no
collisions. Its narrow-view sampling is too noisy for a precise comparison:
decay standard errors reach 0.297/us. The original automatic log-range plot is
poor at displaying nearly empty bins; the inspected `_linear.png` displays the
same saved data without changing a value. Do not treat these noisy fits as a
successful experimental validation or discard inconvenient batches.

The improved calculation integrates initial axial position exactly. Each random
history carries the interval of initial source positions that has not yet left
through an end. Every flight intersects this interval with the finite cylinder;
all lost intervals retain their actual exit-time range. The detector scores
the remaining source measure visible through the slit. This removes rare
slit-hit noise without widening the aperture, discarding photons, assuming an
infinite cylinder or renormalizing energy. See `resonance_axial.md`.

Run 95409 exited 0 after 575.457 s: a 65536-history vacuum trajectory and 32
complete absorbing batches (4096 histories each), in
`outputs/resonance_axial_experiment`. Source:
`1100b61556b8529eb77f57f92e4bf0fa241d4ebaaa0e862bab2ad09399d77e4a`.
Both collision and no-collision cases use 293 K, eight independent seeds per
pressure/model and NIST natural rates. PNG inspected; PDF, full exit intervals,
source/input archives and detector histories saved. No unit suite was used.

## Collision physics and independent limits

The report's printed p35 gives nu_c=0.904 A (gu/gl) N lambda^3/(6 pi^2).
Lorentz frequency HWHM becomes (A+nu_c)/(4 pi). An excited atom still waits
Exp(A); conditional on its actual residence time t, a decorrelating collision
occurs with probability 1-exp(-nu_c t). This preserves the collision/time
correlation rather than independently flipping a fixed probability. A reset
redraws thermal velocity and rest frequency; it destroys frequency memory,
not the excitation. It is not quenching or additional radiative decay.

Rates are 1.3594931e7 and 3.8842660e7 /s at 14.7 and 42 mTorr. The analytic
integrated reset probabilities are 0.02491763 / 0.06804442; measured weighted
batch fractions span 0.02483536–0.02505181 / 0.06785467–0.06809624. These are
observations of the prescribed model, not independent collision cross sections.

The integrated-source vacuum detector integral is 5.684105110424832e-5/sr,
versus slit/(4 pi L)=5.684105110424834e-5/sr. Detector CDF differs from independent
disc quadrature plus natural decay by at most 0.00153536. Independently averaging
ray paths over beam position, direction and uniform axial source gives side
escape fraction 0.85428718 versus 0.85425088 ± 0.00067990, and mean escape
1.92792360 ns versus 1.92659549 ± 0.00732073 ns. All absorbing histories finish;
latest source exit is 67.4741 us. Maximum source-weight budget residual is
1.33e-15, with no source-weight repair. This is constant-hnu0 inventory
conservation, not a coupled plasma/neutral/photon energy theorem.

## Experimental result and remaining discrepancy

Rates below are inverse microseconds; ± is one whole-batch jackknife standard
error. Published bounds are ±3%, not specified standard deviations. The paper
does not give exact numerical early-fit intervals. All three declared windows
are retained; no window or atomic rate was selected to force agreement.

| Pressure | Fit window | No collisions | Resonant collisions | Published early rate |
|---|---|---:|---:|---:|
| 14.7 mTorr | 0.25–2 us | 0.37499 ± 0.01474 | 0.45685 ± 0.01495 | 0.468 |
| 14.7 mTorr | 0.5–3 us | 0.36951 ± 0.01570 | 0.45372 ± 0.01112 | 0.468 |
| 14.7 mTorr | 1–5 us | 0.38803 ± 0.01086 | 0.42860 ± 0.01046 | 0.468 |
| 42 mTorr | 0.25–2 us | 0.20343 ± 0.02274 | 0.23842 ± 0.01405 | 0.310 |
| 42 mTorr | 0.5–3 us | 0.15622 ± 0.01699 | 0.27648 ± 0.00945 | 0.310 |
| 42 mTorr | 1–5 us | 0.14908 ± 0.00977 | 0.26009 ± 0.00946 | 0.310 |

Collisions improve the early comparison, particularly at 14.7 mTorr, but the
42 mTorr result remains 10.8–23.1% below the reported early rate. The component
is not experimentally validated across both pressures. Global mean escape
times are 2.63164/2.40192 us at 14.7 mTorr and 6.11041/4.08850 us at 42 mTorr
(without/with collisions); their inverses are not the detector observable.

Remaining physical uncertainties include actual beam/pulse shape, finite
angular/spectral optical acceptance, quenching, polarization, atom displacement,
and the collision prescription. Formal Voigt tails reach relative detuning
0.11544 in this sample, outside the constant-energy narrow-line approximation.
They remain in the evidence; no frequency cutoff or photon correction is used.
Their spectral contribution still needs quantification. Plot reference curves
use the published slope with an illustrative amplitude; no absolute measured
photon-count amplitude is claimed. Impurity-flagged table rows were not used.

The report also contains the original channel counts for runs151/197 on printed
pages124/127 (PDF135/138), channel width0.3289us. Next inspect the instrument
time convention and those raw curves, then address actual source/optical and
spectral assumptions. Do not tune nu_c or declare agreement from one window.
Keep frc-transport.5.1 and its parent open; moving atomic/radiation coupling
remains absent and must use physical collision/excitation data.
