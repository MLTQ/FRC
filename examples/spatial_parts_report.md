# Spatial contribution and physical charge report

Reads only completed, replayed contribution, physical and one-sided charge
audits. Draws signed physical Btheta directional contributions for the saved
error and independent mode, normalized volume/trace charge defects, and the
known-mode per-term versus full axial-current errors.

Projection bars have units1000/s and can cancel; they are not eigenvalues or
percentages of a cause. Charge diagnostics are dimensionless and distinguish
volume divergence, radial jumps and axial jumps. Different current-error bars
use each component's own continuum target norm, explicitly labeled.

Archives exact input hashes, source/companion files and summary, then writes
PNG/PDF. Visually inspect the PNG. PDF is the same figure, not independently
inspected. No simulation is rerun or changed, and no unit tests are used.
