# Experimental boundary-order probe

A declared fixed strip selects first-order shared face/corner fluxes and metric
sources incident to the wall/end cells, while retaining the second-order
interior. The existing conservative assembly, CT, natural boundary work and
conditional adaptive fallback then apply unchanged. This is a numerical
order reduction with vanishing strip width under refinement, not new material
physics, an atmosphere or a thermal repair. Production defaults are untouched.

`install_probe` wraps only flux selection on one solver instance, recording the
forced and neighboring affected cell counts. Existing limiter diagnostics count
only additional adaptive fallback; zero there must not be called full second
order. Installing twice or on an already first-order solver is rejected.

`run` records unchanged shock inputs, reference, source, complete controls,
history/profiles, accepted failures and an explicit experimental label. Engine
and probe source hashes are separate. Partial files are marked incomplete.
This controlled hypothesis must pass source-driven flow and boundary/mesh/time
studies before adoption; a reduced global energy error alone is insufficient.
Earlier disposable nr192 exploratory runs are not reused as fingerprinted data.
