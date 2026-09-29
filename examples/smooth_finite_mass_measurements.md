# Complete smooth finite-mass measurements

Freezes source before all runs, retains a source archive, then records full inputs,
rate/face/source arrays or initial/final trajectories, parameters, configurations,
accepted histories and exact NPZ hashes. Partitioned runs use disjoint protocol
indices and must share the same frozen source. Initialize the shared source
manifest before parallel writers begin; a source change is an error.

Rates use independent parent species or exact contact pressure/energy targets.
Trajectories retain positivity, conservation/source budgets and exact accepted
times. Failed steps preserve fields, all diagnostics and attempted dt without
retry; they remain partial. Initialization/rate failures are explicitly distinct
from accepted trajectories and retain input state when available.

Contact initialization errors are separately recorded. Short-contact changes
are normalized to the actual analytic change; at a full transit the true change
is zero, so its relative denominator is explicitly null and dimensional error
is retained. Wave perturbation norms and manufactured change norms retain the
baseline independent definitions. Nothing is normalized to hide small signals.

Induction records distinguish strict remainder solve residuals, full-equation
backward errors and descriptive rate-relative conditioning. A successful run is
not an independent audit or a formation claim; a read-only report must verify
all configurations, source/field hashes, physical accounting and analytic norms.
Completion is serialized as a native boolean even when a computed timestep makes
the time scalar a NumPy value. The first attempt's output failure and orphan
contact fields are preserved separately; they are not claimed as audited runs.
