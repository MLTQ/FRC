# Finite-cell detector comparison with exact axial source averaging

`python -m examples.resonance_axial_experiment NEW_OUTPUT --packets 4096 --batches 8`
compares complete finite-cylinder photon histories with clean ORNL-TM-4368
14.7 and 42 mTorr detector rates. It includes both zero-collision controls and
the report's resonant decorrelation/broadening prescription, with NIST atomic
rates. Both use the report's 293 K temperature. No experimental rate is fitted
to set physics. The beam disc remains an assumption based on the collimator.

The transport integrates uniform initial axial position exactly, while retaining
all finite-end source losses and their actual exit times. This removes the rare
slit-hit noise found in the direct-position study. Each seed gives an independent
batch of complete random histories. Detector batch SEM and three-window decay
jackknife uncertainty include within-history correlations. All source weights
exit; no cutoff, weighting repair, or stochastic killing is used.

A complete 65536-history vacuum calculation compares projected detector
integral/CDF with independent disc quadrature and exponential decay. Independent
beam/direction quadrature also predicts side fraction and mean escape time:
for side-ray distance d and |direction_z|=mu, averaging a uniform initial z gives
P(side)=max(1-d*mu/L,0), and E[path]=d-mu*d²/(2L) if d<=L/mu, otherwise L/(2mu).
This reference does not use the transport's evolving surviving-source intervals.

Fit windows, rates and optical limitations retain the conventions documented
in `resonance_experiment.md`. Published detector rate bounds are not standard
deviations. Plotted published exponential amplitude is anchored for illustration,
not used to modify simulated data or extract the simulated rate. Photon/source
budgets are distinct from missing plasma energy, neutral motion and quenching.

Full terminal exit intervals, source-weight budgets, collision statistics,
detector histories, source/input archive and scientific figures are retained.
Disagreement or fit-window dependence is a result, not a reason to tune a rate.
No unit tests or new reporting framework are introduced.

Run95409 completed all33 trajectories in575.457s. See
`docs/RESONANCE_EXPERIMENT.md`: the14.7mTorr early comparison improves with
sourced collisions, while42mTorr still differs10.8–23.1% over the declared
windows. This is not broad experimental validation. Formal tails reach11.5%
relative detuning and remain a stated narrow-line/optical-response limitation.
