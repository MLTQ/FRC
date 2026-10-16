# Argon cell and narrow detector signal

`python -m examples.resonance_experiment NEW_OUTPUT --packets 4096 --batches 8`
follows complete 104.8 nm photon histories in the ORNL-TM-4368 cell and compares
the narrow side-view signal with clean Table 1 rows at 14.7 and 42 mTorr.
Geometry is 7 cm long, radius 1.11 cm, 0.714 cm diameter aperture and 50 um slit.
Initial excitation is uniform along the axis and a 3 mm diameter disc. That
disc is an assumption based on the upstream collimator, not a measured profile.
The pulse is impulsive; finite source timing, atom motion and optical response
are unresolved. Temperature is assumed 300 K. No rates are fitted to the data.

Source z is importance-sampled from 25% uniform plus 75% truncated Gaussian
of sigma 3.7 mm, with the exact likelihood ratio `(1/L)/q(z)` scoring the
detector. There is no source renormalization. The positive uniform mixture
retains support everywhere in the real cell. Every analog packet is followed
to a real side/end exit, including after the 10 us detector plot interval.
Counts in the transport file count computational histories; the weighted
detector signal estimates the physical uniform source. Constant h nu0 energy
and all prior transport exclusions remain explicit.

Independent seeds define whole-history batches. The signal is a sum of next-
event transmission probabilities per steradian, with the actual conditional
observer frequency. It is not a histogram of photons forced toward the detector.
Variance and leave-one-batch-out jackknife uncertainty preserve correlations
between repeated emissions of a packet; emissions are not independent samples.

Three declared fit windows, 0.25–2, 0.5–3 and 1–5 us, remain visible. The rate
matches the conditional mean bin time for equal-width exponential bin
probabilities, with free amplitude. Negative apparent decay is not clipped.
The paper gives a window-selection procedure but no exact windows; matching
one chosen interval is not validation. The plotted reference line has one
amplitude anchored over 0.5–3 us for display only; slopes are not adjusted.
Published ±3% bounds are not relabeled as statistical standard deviations.
Impurity-flagged rows remain separate and are not used as clean comparisons.

A complete 262144-packet vacuum run compares detector integral and the arrival
histogram against independent beam-disc quadrature plus exponential lifetime.
For this beam fitting fully inside the vertical aperture, the exact integral
is `slit_width/(4*pi*length)`. This checks the source/view geometry and likelihood
weighting, not many-scattering experimental accuracy. No unit tests are added.

The current calculation omits pressure broadening and quenching, polarization,
recoil, frequency-dependent optics, and a measured beam/pulse profile. Disagreement
is preserved as a physical-model discrepancy, not corrected with an escape
factor. Results are a component comparison, not a plasma formation prediction.

| Dependent | Contract |
|---|---|
| `resonance_transport` | supplied in-gas source positions; full terminal histories |
| `SideView` | likelihood weights and fixed +x narrow-aperture radiance |
| ORNL archived CSV | measured early detector slope, pressure in Torr, impurity flags |

Run27479 completed all17 trajectories. Direct-source slit sampling was too
noisy; preserve its data and use `resonance_experiment_linear.png` for a readable
display. The subsequent exact axial-source calculation retains finite ends and
reduces this variance. See `docs/RESONANCE_EXPERIMENT.md` for both results.
