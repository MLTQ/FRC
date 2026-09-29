# Current inertia verification

Finite Maxwell/Faraday derivatives check the independent continuum sphere
mode and physical Ohm law. Actual-mesh current/local heat converge against it;
a separate modified-Bessel collisionless sphere checks source screening and
the vanishing-carrier current limit. Neither reference uses solver matrices.

An independently derived variable-density/resistivity sphere checks continuum
current and local-heat convergence with correlated, explicitly prescribed
coefficients; its mathematical construction is not a species closure.
Independent physical-volume magnetic/current energy checks nonuniform density,
resistance and deformed geometry, source work and per-cell positive heating.
The explicit massless limit matches the old resistive component. A separately
diagonalized spatial mode isolates midpoint time order. Invalid density and
failed time proposals do not change caller fields; exact zero state remains zero.
The exactly ideal massless limit preserves arbitrarily small material flux
under a much larger bias change; no roundoff-generated flux replaces it.
Its independently reported Ampere constraint is checked as well.
Stationary nonzero collisionless current is unchanged exactly by the increment
solve and passes the backward-residual check. Current convergence retains the measured initial weak-derivative surface
error (approximately order 1.5), which inertia no longer rapidly damps. The
massless heat comparison uses a global heat-scaled roundoff allowance for cells
whose true heat is tiny, while retaining a tight relative check elsewhere.
