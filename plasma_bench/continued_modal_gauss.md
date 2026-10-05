# Continuation proposals for the same Gauss equations

Adds two optional starting/correction proposals to the frozen accelerated
Gauss algorithm. The previous collocation polynomial is extrapolated to the
next two stage times; accepted residual secants may be reused. Both are used
only when dt and the supplied physical starting state/field exactly match the
previous successful endpoint. A changed input or time step starts fresh.

For Gauss nodes c1,c2, integrated Lagrange weights b1(s)=(s²/2-c2 s)/(c1-c2),
b2(s)=s-b1(s), the next stage change proposal is
h sum_j [bj(1+ci)-bj(1)] Fj(previous). This is a guess, never a solution or
an external forcing. The actual rounded full nonlinear residual is evaluated
before acceptance. All original chord/Newton fallbacks, positivity, induction,
endpoint and full-stage bounds remain unchanged. All Fourier modes are kept.

Continuation memory commits only after a successful endpoint; local list
copies prevent a failed trial from changing previous history or input state.
The `continued` diagnostic records use. Stage histories reveal the proposed
initial residual. The explicit step remains isolated from previous frozen
source closures; it adds no physical damping, projection or energy repair.
