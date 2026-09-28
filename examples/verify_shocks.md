# Shock-verification report and figure

Builds the full measurement suite, comparison table, six-panel scientific figure
and Markdown report. Retains the centered method's false reversal and irregular
frozen-ratio convergence, separate numerical heating, time/coefficient controls,
boundary contamination checks and direct/Picard agreement. Local profiles are
compared at identical radial points; controls state their shorter comparison
time explicitly. Data and report-source fingerprints are separate so formatting
can be rerendered without rerunning unchanged expensive physics measurements.
Boundary sensitivity is measured at nr=48; the report does not extrapolate that
control into a guarantee for every grid or longer evolution.

`python -m examples.verify_shocks --reuse` loads only per-case files whose
measurement/package source, environment and fixed case specification match.
Changed physics invalidates reuse. Failed cases remain visible and include
accepted-state diagnostics; they are never presented as completed runs.
