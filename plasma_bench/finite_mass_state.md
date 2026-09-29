# Moving finite-mass state and thermal partition

Stores five conservative primal cell variables: rho=(mi+me)n, radial/axial
momentum, total energy and axial angular momentum; corner psi, primal g and
dual electron U are separate. Both species gammas are explicit. No old MHD
state or energy convention is silently reinterpreted.

Evaluation maps the actual total mass to dual n, computes compatible magnetic
and physical reduced-mass current energy, electron energy from Q1 pe, and
recovers ion energy by subtracting these and bulk kinetic energy from total E.
This subtraction defines the partition; it never changes total E. Ion pressure
is overlapped onto the dual grid for Ohm/enthalpy terms. Distinct primal and
dual pressures and their interpolation error must be verified.

Bulk radial/axial velocity uses cell momenta divided by rho; rotation uses
L/(rho <r²>), with toroidal kinetic energy rho <r²> Omega²/2 as in ToroidalState.
Radial subcell kinetic variation is not recovered from its cell mean; this is
an explicit spatial approximation. from_primitive initializes the numerical
energy convention once; an independent continuum initializer/norm must measure
its representation error. No accepted step resets energy to these components.

All density/species states must remain strictly positive; nonfinite/negative
states reject without floors. This is a filled-cylinder reference state, not
an unconstrained magnetic domain: regular-axis, fixed-wall and periodic flux
conditions are checked during evaluation. The model is not
an exterior vacuum carrier prescription. Parameters support physical electron
mass or an explicit consistent massless comparison, no effective mass fitting.
