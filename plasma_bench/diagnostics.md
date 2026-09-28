# Physical diagnostics
`characteristic_scales` evaluates explicitly labeled reference scales in SI
from supplied material parameters and |bias|. It returns finite JSON-compatible
numbers or null for undefined quantities and specific missing-physics warnings.
Drude collision/Hall estimates are inferred from arbitrary sigma and are not
physical transport predictions. `magnetic_energy` integrates total field in a
specified region; it is distinct from the solver's exact induced-energy budget.
`longest_sampled_interval` reports only the lower bound between consecutive
positive topology samples, never inflating an isolated positive frame's lifetime.
Waveform feature times, quasistatic light-transit ratios and oscillatory skin
depths trigger explicit resolution warnings with roundoff-tolerant thresholds.
Sampled-wave knots count as features.
