# Polarized resonance transport against actual channel counts

`python -m examples.resonance_channel_experiment NEW_OUTPUT --packets 4096 --batches 8`
runs the finite ORNL cell with the sourced collision rate, comparing isotropic
emission with polarization-preserving electric-dipole scattering. Geometry,
atomic rates, neutral temperature, source distribution and collision frequency
are identical. Both calculations finish every source history. No rate or
collision parameter is fitted to obtain agreement with the experiment.

The unaligned initial dipole and complete orientation reset at a resonant
collision are explicit assumptions. This is a zero-field isolated J=0-1-0
limit; initial proton-excitation alignment, depolarizing cross sections,
instrument polarization sensitivity and quantum line interference are unknown.

Raw counts and metadata are visually transcribed in data/atomic/radiation_reference.
Reverse timing is retained. The first full decay bin starts an unknown offset
after excitation; offsets0,0.5,1 channel are all compared. Detector scores are
accumulated on exact half-channel edges and folded by the32us repetition period.
Actual photon histories never wrap. The pulse width (reported <=20ns FWHM) and
electronics response are not known per run and are not deconvolved or fitted.

Fixed windows contain6,9,15,24 complete bins (1.9734,2.9601,4.9335,7.8936us).
Measured counts use a Poisson likelihood with the independently reported
background held fixed. The exponential's amplitude and slope are fitted only
to summarize each measured window; the conditional expected-Fisher slope error
does not include background uncertainty or model departure from an exponential.
Simulation slopes use their own detector shape; whole-batch leave-one-out
jackknife errors retain within-history correlations. All offsets are preserved.

For direct count overlays, only one positive multiplicative intensity is fitted
to the simulated shape, accounting for unknown excitation yield/detector
efficiency. No time shift, shape, atomic lifetime or collision rate is fitted.
The reported Poisson deviance is conditional on the Monte Carlo mean, not a
calibrated goodness-of-fit p-value. Display bands show batch SEM at that fitted
intensity, excluding its uncertainty. Observed spectral-tail contributions are
retained without imposing an unmeasured monochromator passband.

The same end-to-end run includes vacuum transport against independent source/
direction quadrature and natural decay, plus analytic dipole angular moments.
The scalar and dipole trajectories retain complete exits, Stokes scores, source
budgets and input/source identities. This is a concrete physical extension and
experimental comparison, not a test framework or moving-plasma validation.

Run67602 completed33 trajectories in611.819s; results and inspected figure are
in outputs/resonance_channel_experiment. Dipole scattering does not resolve the
42mTorr mismatch. A rare formal relative detuning1.00811 exposes the inherited
narrow-line approximation explicitly; no clipping was used. See the new top of
docs/RESONANCE_EXPERIMENT.md for all windows, uncertainty and source budgets.
