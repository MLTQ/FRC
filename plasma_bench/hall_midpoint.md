# Nonlinear Hall/resistive implicit midpoint

Holds carrier density and fluid momentum fixed for this induction substep.
For midpoint current j, electric coefficients e, and magnetic Hessian H:

    M j + (dt/2) C.T H C e = -C.T H x_old
    M e = (X(B_mid)+R) j,     x_new=x_old+dt C e.

Each Picard iteration updates B_mid and solves this sparse mixed system;
the stiff linear current response is implicit. A physical impedance scales
e into current units, and diagonal mass scaling balances cylindrical bases.
An existing sparse LU is only a preconditioner for the current, unchanged
linear system; GMRES checks the actual residual. If it fails its bounded
iteration budget, a fresh exact-system LU is formed. Reuse changes cost, not
the Hall equation. This computational cache may change on rejected trials.
The final nonlinear Ohm residual and linear solve residual must pass explicit
tolerances. Iteration failure raises without mutating input fields. No damping,
resistivity adjustment or state repair is introduced to obtain convergence.

At convergence the quadratic energy change equals minus dt j.T R j, while
j.T X j=0. Positive physical eta |J|² is diagnosed once. Total fluid energy
receives only midpoint Poynting divergence; local transfer errors are retained.
The result carries new psi/g, that conservative energy change, physical heat,
iteration/residual counts and electric coefficients. Shapes and finite initial
axis/wall/periodic flux are validated before solving or taking a zero-transport
shortcut; the Faraday map preserves those traces. This does not by itself
advance mass or momentum.

Midpoint avoids an explicit whistler stability restriction but still has phase
error and nonlinear convergence limits. A large energy-conserving step is not
accurate wave propagation; mesh/time/solver-tolerance controls are mandatory.
