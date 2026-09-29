# Physical trajectory comparisons

Directly evaluates each mesh's physical B, current basis and Q1 species pressure
at Gauss points on fine-cell subdivisions. Common-time nested meshes and fixed-
mesh time/tolerance/quadrature controls therefore compare actual represented
fields, not restricted coefficient arrays alone. RMS differences retain SI
units. Relative differences use the finer trajectory's actual change and its
initial field separately; unresolved zero change gets null, not a false norm.

Time ratios use a stated 320-step numerical reference, not an exact continuum
claim. Nonlinear mesh differences are also inter-trajectory measurements, with
no exact nonlinear solution implied. The independently exact wave remains the
continuum spatial check. Comparison requires completed equal-time trajectories
and cannot quietly include a failed partial run.
