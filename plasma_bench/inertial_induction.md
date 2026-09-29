# General instantaneous finite-inertia induction

Solves the complete poloidal/toroidal rate system using the existing compatible
electric mass M, physical inertia I, magnetic energy H and Faraday C. Fixed
axis/wall flux degrees have zero rates and are eliminated explicitly. On the
remaining magnetic space H is positive definite. The equivalent SPD equation is

    [H + H C M^-1 I M^-1 C.T H] Bdot = H C M^-1 f.

Diagonal scaling balances flux/toroidal units. The exact scaled-H factor is a
preconditioner only; conjugate gradients solve the full unchanged operator.
Reconstructed Ohm, Ampere, Faraday and operator residuals must pass. Outputs are
checked in both the original scaled dual residual and the physical magnetic
energy norm for Faraday. Up to three explicitly reported linear iterative
refinements solve the original reconstructed Faraday residual. These are
corrections to one unchanged linear problem, not physical time-step retries;
the requested tolerance is never relaxed. Iteration totals and each residual
component are returned. Failure after that budget remains atomic.
Outputs are rates; input state/forcing is never changed. Model density/inertia remain fixed
for this instantaneous solve, but forcing may include the actual n_t required
by a moving system. This is not a frozen-density time-step approximation.

InertialInduction is tied to its model/mesh and positive carrier domain. It
supports an explicit consistent massless model. It neither adds vacuum carriers
nor modifies physical parameters when convergence fails. Direct mixed-system
comparisons and independent analytic/manufactured fields check the reduction.
