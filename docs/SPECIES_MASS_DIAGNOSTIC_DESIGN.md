# Entropy-rate mass and representation diagnostic

Completed in this continuation: SPECIES_INTERFACE_VERIFICATION.md records the
full replayed mass split and subsequent F8/F16/F32 fixed-current refinement.
The exact-mass diagnostic alone leaves most of the source unresolved. The next
full-operator control is COUPLED_FLUID_REFINEMENT_DESIGN.md.

The complete q12/16/20 weak interface budget now reproduces original K to about
1e-12 relative. Its production entropy rate differs from the volume-advection
covector divided by the actual fluid volume by only about 2e-8 in the saved
error direction, while the latter still differs about 55% from the strong
physical source. Large face and transfer terms cancel. No single omitted
interface term or entropy-quotient difference explains that discrepancy.

Next compare the actual diagonal point volume V with the independently
integrated physical mass M=S^T W S. For g=-a/e J dot grad(S s0), retain:

    c_lumped = V^-1 S^T W g,
    c_exact  = M^-1 S^T W g.

These are diagnostic vectors only. The original production Cdot is never
replaced or evolved. Reconstruct f=S c/(S n0) and measure the actual field
identity

    f_production - g/n0 = (f_production - f_lumped)
                       + (f_lumped - f_exact)
                       + (f_exact - g/n0).

Save all vectors, physical component fields or fingerprints, pairwise signed
integrals and their sum. Do not infer cancellation from norms alone. Retain the
nonrepresentable strong-source remainder separately from the mass-map defect.
Compare both species and directions at independent quadrature 16/20 and fully
replay. Measure native and physical entropy integrals with a dimensionally
matched absolute-integral normalization; large arbitrary direction normalization
can otherwise make floating residuals look like physical charge/entropy losses.

Independently verify the radial and axial mass matrices against analytic basis
integrals. In each degree-p radial element, lower Legendre degrees integrate
identically under Gauss–Lobatto weights and the physical x integral. For degree p,
the exact squared integral is 2/(2p+1), while the point-weight integral is 2/p.
Their ratio is p/(2p+1)=5/11 here. The real Fourier Nyquist cosine has continuous
mean square 1/2 versus nodal mean square 1, giving the axial ratio 1/2. Compare
the complete normalized matrices with these orthogonal projectors and measure
the actual source content in the highest radial and axial subspaces.

This can identify the accuracy contribution of the declared discrete metric;
it does not prove a cause of coupled rapid growth or authorize swapping one
mass matrix in an otherwise coupled thermodynamic/energy system. Any subsequent
change requires rederiving all paired work identities and complete trajectory,
convergence and conservation checks. Original tolerances/masses/states remain.
