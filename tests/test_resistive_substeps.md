# Stiff magnetic-decay accuracy and atomicity

An independently diagonalized positive lowest spatial eigenmode has an exact
semidiscrete decay exp(-4) over the requested interval. A single raw midpoint
step reverses its sign despite a valid energy/heat balance. Adaptive field,
current and local-heat checks must instead give positive amplitude, measured
error reduction with tolerance and unchanged combined energy without repair.

This isolates time integration from spatial error; continuum accuracy remains
the separate sphere-decay test. Finite tolerance cannot guarantee the sign of
arbitrarily tiny fields or a topology theorem. Exhausting the trial budget
raises without mutating caller data; full-plasma rollback is checked elsewhere.
Finishing on the final permitted trial must succeed rather than falsely exhaust
the budget.
