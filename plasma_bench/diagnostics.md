# Physical diagnostics
`characteristic_scales` evaluates explicitly labeled reference scales in SI
from supplied material parameters and |bias|. It returns finite JSON-compatible
numbers or null for undefined quantities and specific missing-physics warnings.
Drude collision/Hall estimates are inferred from arbitrary sigma and are not
physical transport predictions. A separate `coulomb_reference` now records
conditional Maxwellian collision/mean-free-path, electron-inertia and carrier
scales from supplied temperatures/species, with sourced scope and review flags.
It uses plasma radius, shortest waveform scale and zero initial current; it is
not a substitute for evolving local gradients/current. Original total physical
collision/mean-free-path keys remain null because no complete transport model
has been selected; prescribed conductivity is never overwritten.
`magnetic_energy` integrates total field in a
specified region; it is distinct from the solver's exact induced-energy budget.
`longest_sampled_interval` reports only the lower bound between consecutive
positive topology samples, never inflating an isolated positive frame's lifetime.
Waveform feature times, quasistatic light-transit ratios and oscillatory skin
depths trigger explicit resolution warnings with roundoff-tolerant thresholds.
Sampled-wave knots count as features.
