# open_ends_verification.py

`python -m examples.open_ends_verification` runs a uniform magnetized axial flow
(ρ = 1, p = 1, v_z = 0.5, uniform field) with open ends at 16 and 32 cells,
plus 16 cells with dual energy. Downstream of the upstream rarefaction head
(at z = −1 + (v + c_f)t) the state must stay exactly uniform. It reports that
error, the outflow mass against the downstream exact ρ v πR² t, and the mass
and energy residuals.

`--conduction` checks the open-line conduction loss against an independent
ODE integration in a uniform resting plasma.
