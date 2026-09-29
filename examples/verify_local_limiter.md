# Local limiter report and scientific figure

Consumes the independently persisted local/global study. Shows original shock
error/profile changes, the actual sparse limiter support, signed independent
response-energy errors, and source-work refinement. Failed runs are excluded
from final-time curves and recorded explicitly with accepted times. Controls
use active-limiter grids; historical inactive cases retain their fingerprint.

The report distinguishes exact fluid flux conservation from full-space energy
accuracy and explicitly reports the latter's nonmonotonic shock error. It does
not interpret positive pressure, sharper shock profiles or a passing test suite
as validated formation. Renderer source has its own fingerprint, so report
changes do not invalidate the numerical measurements.
Shock controls report both profile changes and independent energy residuals;
pulse controls separately quantify changes in work and peak speed.
The report also generates complete, source-matched instantaneous energy audits
at nr=96/192, available separately with `python -m examples.audit_interface_energy`. Their
independent perturbation check validates the rate identity while exposing
spatial error; it is never presented as a solver-accuracy or energy-repair test.
The shock residual is also compared with net outward fluid energy, explicitly
distinct from source work, to avoid hiding interface error behind initial energy.
Plotting alone can preview available completed cases for layout inspection;
the final `verify` path still requires every configured case and both audits.
With `--reuse`, missing, incomplete or source-mismatched audits are recomputed;
a fresh checkout can reproduce the whole report with one command.
The nr=384 outcome is discussed explicitly: interior error refines while the
independent interface energy error and numerical heat worsen. Its runtime,
steps and minimum density are recorded; no density floor hides the thin layer.
Large scalar-array trees use lossless compact JSON formatting, checked by a
round-trip parse. Object records remain readable and the formatter has its own
fingerprint, separate from numerical provenance.
