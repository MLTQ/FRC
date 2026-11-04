# HOPS 25 µs release configurations

`hops_release_25us_pulse.json` commands one +3 kV, 200 ns pulse and then release.
`hops_release_25us_zero_bias.json` is the matched zero-voltage control; only bias
and label differ between them. Both retain the same ideal electron emitter
until 200 ns and stop it afterward. No subsequent emitter or electrode pulse.

These files are complete config inputs to `examples.iec_3d`. Their 40 µs fixed
period and 0.005 duty are an implementation of a single 200 ns pulse over the
25 µs interval, not a proposed resonant frequency. A longer duration beyond
40 µs requires a new pulse protocol; never simply extend the duration.

Both use 5 ps integration, 25,000 saved intervals (1 ns cadence), unchanged
18/1 mm mesh and original seed/plasma/source inputs. The fixed mode ignores the
legacy modulation/frequency fields retained for config compatibility.

The complete physical protocol, interpretation and execution limits are in
[HOPS_RELAXATION_RUN.md](../docs/HOPS_RELAXATION_RUN.md). Cases are designed,
not launched. Keep source-frozen provenance and original runs; a full-state
prefix comparison is required during execution. No existing checkpoint may
be resumed with these changed configs through the strict current loader.
