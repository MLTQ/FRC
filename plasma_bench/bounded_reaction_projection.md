# Conservative bounded reaction moments

`bounded_moments` supplies an explicit alternative transfer for the precise
moving-material solver. Fine cell numbers, physical masses, supports and all
inherited velocity knots are unchanged. It remaps momentum and total energy
together; it never changes a source population or clips a thermal state.

For fine cell c, let m_c be its scaled physical mass, P_c/E_c its quadrature
moments and a_ci the cell-average nonnegative P1 shape. The low-order nodal
states have mu_i=sum_c m_c a_ci, P_i=sum_c a_ci P_c and
E_i=sum_c a_ci E_c. They are convex averages of admissible cell donor states.
They conserve total momentum and energy; mass remains the existing fine P0
representation. Bounds on each velocity component and specific total energy
are the incident donor-cell extrema. These cell-mean bounds are deliberately
diffusive at smooth extrema and boundaries; affine exactness is NOT claimed.

The high-order target solves the existing consistent matrix for velocity and
specific total energy. Antisymmetric adjacent-node fluxes consist of
M_ij (z_i^H-z_j^H) plus the quadrature-minus-cell-mean load. A single correction
factor per edge applies to all momentum components and energy. Linear bounds
and the convex constraint e >= |v|²/2 restrict this factor; thermal feasibility
uses 54 bisection iterations. Edges are processed left-to-right, once. This
ordering can introduce asymmetry and must be measured, not called invariant.
These are conservative flux corrections, not state clipping or a global reset.

Return to fine cells uses the same positive shares. With nodal specific heat
theta_i=e_i-|v_i|²/2 and exact cell-average shape product <phi_L phi_R>,

    U_c=m_c [a_cL theta_L+a_cR theta_R
             + (1/2)<phi_L phi_R>|v_R-v_L|²].

This expression is nonnegative in exact arithmetic. Adding the actual
consistent P1 kinetic energy gives sum_i mu_i e_i globally. The P1 velocity
has momentum sum_i mu_i v_i. Local total energy is redistributed; all thermal
changes and redistribution are numerical transfer effects, not atomic heating.
Roundoff causing negative recovered heat rejects; no floor is applied.

The method is motivated by conservative bounded-remap literature, particularly
[Kuzmin et al. (2010)](https://doi.org/10.1016/j.jcp.2010.08.009), but the above
cell/node distribution and ordered correction are explicitly this implementation,
not a claim to reproduce that paper's synchronized failsafe algorithm.

Contract: the precise projection passes cell-endpoint integrated momentum and
energy loads in its per-species normalized coordinates. Returns nodal velocities,
fine thermal amounts, signed fine total-energy redistribution, and correction
diagnostics. Original Galerkin/local-total transfer remains the default.
Full trajectories and independent spatial/source-time comparisons are required
before adoption. Positivity and boundedness alone are not local accuracy.
