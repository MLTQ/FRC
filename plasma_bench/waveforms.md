# Waveforms
`Waveform` supplies vectorized current and time derivative in SI. DC, Gaussian,
cosine-ramped Unipolar, causal DampedSine, Hann-windowed SineBurst and Sampled are
frozen serializable dataclasses. DampedSine has a right derivative at onset.
Sampled clamps endpoint values outside the table, uses right slopes at knots,
and zero derivative outside. No derivative impulse exists because values are
continuous. `Callback` requires both functions and is script-only; serialization
rejects it so metadata cannot falsely promise reproduction. `to_dict/from_dict`
form the explicit JSON schema; no arbitrary code evaluation occurs.
