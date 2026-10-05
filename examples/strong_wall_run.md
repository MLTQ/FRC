# Immutable wall-conforming trajectory runner

Freezes the complete explicit strong-wall trajectory source/runtime closure,
refuses overwrite and preserves copied accepted progress states with raw hashes.
Success saves all final metrics/arrays; a raised failure writes its traceback
and last accepted progress. Uses the separate StrongWallEvolution protocol,
with CLI kind/elements/steps/amplitude; the only time method is full Gauss.
This is candidate verification, not a physical formation recommendation.
