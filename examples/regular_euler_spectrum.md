# Original cylindrical acoustic spectrum in regular variables

Linear equations in dimensionless time ct/R are delta_t=−2∂x(x a),
a_t=−R²∂x delta, where a=R u/c and x=r²/2. The axis a is the dependent
trace of the positive first-element nodes. The matrix is assembled directly
from these equations and shared centered wall/interface fluxes, not obtained
by differencing the nonlinear driver. Physical perturbation energy metric is
diag(Wdx,2x Wdx/R²) with only active a degrees of freedom; it is positive.

Analytic frequencies are positive zeros of J1. The full matrix, metric and
eigenvalues are returned so separate nonlinear RHS and time controls can be
made. Imaginary linear spectrum does not establish nonlinear plasma stability.
Three amplitude controls compare the actual nonlinear driver with this separate
linearization. Real and imaginary mode parts are combined so density **and**
velocity are perturbed, exercising pressure and continuity simultaneously.
Matrix exponential and RK4 half/quarter steps cover two periods.
Frequencies are the sorted positive spectrum, without discarding an inconvenient
branch: degree3's third sorted mode is currently not the analytic third Bessel
mode, and that discrepancy must be retained rather than called verified.
