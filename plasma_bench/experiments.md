# Experiments and parameter sweeps
`comparison_config` normalizes a unipolar, damped sinusoid or Hann-windowed burst
to an equal peak absolute vacuum perturbation at the chamber center over the
run window. Dense bracketing plus scalar optimization finds the peak; source
amplitude is linear. It does not equalize energy or positive excursion.
`compare` saves three complete results. `sweep` varies conductivity and peak/bias
ratio and records complete configurations, warnings, classifications and code
provenance; other parameters can be varied with dataclasses.replace. Zero-center
coil combinations and zero-bias sweep normalization fail explicitly.
