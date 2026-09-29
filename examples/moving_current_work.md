# Current-energy work with changing density

Checks d(WB+KJ)/dt = -j dot f - 1/2 integral I_local |J|² n_t/n.
The rate includes both j dot I J_t and the actual time derivative of 1/n.
The residual is normalized by its component powers, never a fitted heat term.
The mixed-equation residual first verifies the fixed-n part; physical derivative
replay then includes the separately recorded density contribution on both sides.

Separately differentiate physically interpolated magnetic/current energy at two
small intervals using fourth-order Richardson extrapolation of symmetric energy
differences. Magnetic polarization and a common rational denominator for J²/n
evaluate that exact difference before quadrature, avoiding loss of digits from
subtracting background-dominated global totals. Both B and n change
along the computed rate; density variation is never frozen in this replay.
Independent physical quadrature avoids production energy/curl matrices. These
are tangent-energy checks, not proof of a conservative finite time step or a
completed bulk momentum/thermal integration.
