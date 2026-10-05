# Coupled regular moving plasma

## Purpose
Assembles the independently verified bounded components on the actual
n,q,l,mz,Ce,Ci state. This is an unforced filled periodic cylinder, with fixed
axis/wall poloidal flux and impermeable radial fluid boundaries. It is not a
vacuum, coil, end-loss or formation solver.

## Components and contracts
`RegularPlasma.derivative(state, HarmonicMagnetic)` returns a `RegularPlasmaRate`:
fluid and packed magnetic rates, compatible current rate, full electric field,
strict induction residual and grid energy diagnostics. Each evaluation builds
its physical current model from that state's density. Bulk continuity supplies
the single density rate used by current inertia and current advection; Lorentz,
drift stress and species pairs do not change continuity.

The default geometry cache shares only density-independent differential/trace
maps, L2 transfer and magnetic-Hessian factorization. Shallow operator copies
are rebound to the fresh density model at every call; no inertia coefficient,
forcing, solution or density is reused. `cache_geometry=False` rebuilds the same
operators for bitwise changing-density control. Instances are intended for
sequential evaluation, not shared concurrent mutation.

The forcing sums ideal/Hall, changing-density inertia, current advection,
stretching, self drift and reduced species pressure. The compatible species
potential is restored only to full E and the original Ohm residual, not
numerically curled into magnetic evolution. The physical masses and separate
local species temperatures are unchanged. Grid fluid energy plus quadrature
current/magnetic energy is distinct from continuous reconstructed total energy.

No forcing, fitted equilibrium, artificial diffusion, filter, atmosphere,
positivity floor or energy correction is applied. Failure leaves inputs intact.
This module includes no time method or automatic time-step recommendation.
Full coupled trajectory stability, original-equation convergence and physical
interface budgets require independent verification.
