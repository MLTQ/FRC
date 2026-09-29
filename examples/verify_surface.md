# Completed surface study assembly

Requires all 13 independent references and 21 moving cases to be complete and
source/configuration/field-hash matched. Compares actual states on common initial
physical material coordinates, with explicit omitted-volume and boundary-speed
diagnostics. Intersects initial material cells before integrating comparisons,
checks covered physical volume and rechecks every comparison at higher positive
triangle-rule order; diagnostic integration uncertainty remains visible.
Integrations are reused only through a complete exact-input/source-checked cache;
formatting-only edits never require rerunning the plasma or comparison integrals.

Restores final states for diagnostic force decomposition and checks the three
uniform runs against the historical archived arrays, retaining old fingerprints.
Analysis/render sources and Matplotlib version are separately recorded. Writes
the final completed metrics only after both figures and the report are created.
Visual inspection remains required. Initialization in an FRC, positive heat,
smaller cells or a small energy residual cannot unlock formation claims.
The independent initial surface-compatibility audit is separately recorded and
checks its analytic current with one-sided finite Maxwell derivatives. Its
conditional smooth-trace assumptions and primary sources remain explicit.
The final time-control interpretation retains the non-asymptotic velocity ratio
and nonmonotonic current differences; it does not extrapolate a resolved state.
