# Physical trajectory comparisons

Requires matching domain, endpoint and nested meshes. Continuous Q1 electron
pressure differences use exact polynomial integration on fine subdivisions,
in Pa and relative to the finer trajectory's pressure change. Compatible Q1
psi and P0 g prolong exactly to the fine representation; independent physical
B integration measures magnetic differences relative to the fine field change.
These are trajectory differences, not errors against an exact nonlinear time
solution. Representation differences are included in spatial magnetic norms.

Time controls compare pressure 10/20/40 with 320 steps, fast wave 8/16/32 with
256, nonlinear 10/20/40 with 160. Reports both pressure and magnetic ratios.
Opposite wave direction, three finite amplitudes, near-cold electrons and two
half-step standalone diffusion are retained without modifying base cases.
