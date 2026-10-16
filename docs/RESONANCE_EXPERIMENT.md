# Argon radiation: finite-cell detector comparison

## Polarization and original channel counts (2026-10-16)

The radiation solver now has opt-in electric-dipole polarization memory for an
isolated zero-field J=0-1-0 transition. Initial excitation is unaligned; each
resonant collision resets orientation as well as frequency. That orientation
reset is an assumption, not a measured depolarization rate. The scalar model is
retained. Actual and observer rays share the same conditional absorbing atom;
the detector records I,Q,U and spectral-tail contributions. Neither rates nor
photon inventories are fitted or repaired. See resonance_polarization.md and
examples/resonance_channel_experiment.md for equations and limits.

ORNL raw runs151/197 are now visually transcribed in channel_counts.csv with
their original metadata. Reverse timing is documented on printed21-23. Higher
channel is earlier; the first full decay bins are81/77. Run151 channel81 is452,
not the OCR4523 (the last mark is a bracket). Expected background counts per
channel are58.8319875/182.950625 from the report's acquisition/rate/period data.
Absolute time zero was not calibrated. All offsets0,0.5,1 channel are retained;
none was fitted to favor a model. Exact half-channel detector bins are folded
by the32us repetition period; actual histories and exits remain unfolded.

Run67602 exited0 after611.819s:33 complete trajectories, including65536 vacuum
histories and2 pressures x2 scattering laws x8 batches of4096 histories.
All raw exits, detector/Stokes/spectral scores, inputs and frozen source remain
in `outputs/resonance_channel_experiment`; the scientific PNG was inspected.
Source `c4c7fbbdabbd4300e915d9a900838d484606557275b54a3c80072910b43a594b`.
This source change does not enter either live moving-plasma calculation.

Rates below are /us. Data errors are conditional Poisson expected-Fisher
standard errors with fixed reported background; simulated errors are whole-batch
jackknife standard errors. These are not the report's Table1 error estimates.
The simulated values use the nominal half-channel offset; every offset/window
is stored. Longer windows need not share the reported early-time slope.

| Pressure | Full fit window | Raw-count rate | Scalar scattering | Dipole scattering |
|---|---:|---:|---:|---:|
|14.7mTorr|1.9734us|0.51948 +/-0.05189|0.42941 +/-0.01979|0.41697 +/-0.01190|
|14.7mTorr|2.9601us|0.48129 +/-0.03251|0.42813 +/-0.01491|0.42325 +/-0.01143|
|14.7mTorr|4.9335us|0.49452 +/-0.02147|0.42668 +/-0.00577|0.42338 +/-0.00888|
|14.7mTorr|7.8936us|0.49234 +/-0.01796|0.44540 +/-0.00339|0.43731 +/-0.00630|
|42mTorr|1.9734us|0.26705 +/-0.03181|0.27324 +/-0.01920|0.27636 +/-0.02846|
|42mTorr|2.9601us|0.30797 +/-0.01917|0.26768 +/-0.00851|0.27701 +/-0.01526|
|42mTorr|4.9335us|0.33092 +/-0.01128|0.26907 +/-0.00689|0.25559 +/-0.00682|
|42mTorr|7.8936us|0.32050 +/-0.00791|0.26440 +/-0.00283|0.25598 +/-0.00496|

For the42mTorr longest window, all three offsets give scalar0.26333-0.26535
and dipole0.25598-0.26011. Polarization and time-zero uncertainty do not remove
the mismatch. In the count overlays only one intensity multiplier is fitted
(unknown excitation yield/detection efficiency), never a physical rate or time
shift. Conditional Poisson deviances for24 bins are29.79/31.09 at14.7mTorr
and73.71/103.93 at42mTorr (scalar/dipole). These do not include Monte Carlo
uncertainty and are not calibrated goodness-of-fit p-values. Keep both models;
neither is an experimentally validated radiation closure for moving plasma.

Independent dipole direction moments are(0.200219,0.399408,0.400373), against
(1/5,2/5,2/5), with SEM(0.000591,0.000859,0.000859). Unpolarized beam
<cos^2(theta)>=0.398317 +/-0.000860 against2/5. Vacuum detector intensity is
(5.67132 +/-0.00993)e-5/sr against5.684105e-5/sr; Q,U agree with zero within
their sampled errors. Side escape0.854066 +/-0.000677 vs0.854287. Mean vacuum
exit1.909965 +/-0.007269ns differs from1.927924ns by-2.47 standard errors;
retain this fluctuation. Its sampled excited hold is1.861735 +/-0.007268ns
vs1.879699ns; flight time0.0482302 +/-0.0000692ns agrees with0.0482243ns.
The vacuum case was not repeated to select a closer result. Maximum source
inventory residual across absorbing trajectories is1.11e-15.

Dipole time-integrated Q/I is-2.98%/-3.84% at14.7/42mTorr (Q positive along
the tube axis); instrument polarization sensitivity remains unmeasured.
For42mTorr dipole histories, detector intensity fractions beyond relative
frequency detuning1e-4,1e-3,1e-2 are1.628%,0.1887%,0.01685%. These are
sampled integrated contributions, not bounds on every time bin or a measured
spectral response. Formal detuning reaches1.00811 times the line frequency!
That explicitly exceeds the narrow-line, constant-hnu0 approximation; the
saved maximum is unsigned. No frequency or photon was clipped. Small detected
tail weight does not validate a coupled photon-energy budget or the unbounded
emission model. A physically justified finite-frequency kernel and energy
exchange are required before coupling this component to plasma energy.

Next radiation work: resolve source/optical and redistribution assumptions with
independent physical inputs, including the far-wing validity problem; do not
tune rates to the raw counts. Collisional excitation data and moving species/
energy coupling remain open. Pending LXCat terms have not been accepted.

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
