# Independent wall-trajectory physical audit

Uses the established independent root-polynomial quadrature16/20 reconstruction
of original physical density, velocity, both pressure/temperature fields, J/B/E
and conserved energy. Reconstructs StrongWallPlasma only for full E and strict
declared induction diagnostics; the original raw-curl defect and magnetic wall
trace are reported separately. Grid energy and continuous energy remain distinct.
Stable thermal difference evaluation changes no physical state or energy.

Verifies original source/runtime/raw hashes, freezes its own complete closure,
and optionally replays every trajectory array and non-runtime metric bitwise.
No terminal success is written on failed replay. Preserves particle/entropy
changes, point/sample positivity, physical temperature extrema and axis/wall/max
norms. The explicit wall boundary approximation needs mesh convergence; global
energy is not a certificate of local enthalpy/charge interface balances.
