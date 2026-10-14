# Magnetic exterior around the side and both finite ends

FiniteCylinderVacuum statically eliminates a genuinely magnetic-only exterior
around a finite cylindrical plasma. Boundary flux unknowns are the side wall
and both end-cap nodes, excluding the regular axis and duplicate corners.
There are no exterior fluid, conductivity or evolving current variables.
The reduced positive Hessian maps these actual boundary values to vacuum energy
and its work covector. It must be coupled to the same plasma polynomial trace.

Five tensor-product vacuum patches surround the cylinder. Their radial/axial
polynomials match the plasma's traces on each interface; axis flux is zero.
An outer rectangle encloses the whole plasma sphere. Outside it, the existing
decaying multipoles and their all-space energy are retained. Source coils are
separate prescribed fields; these multipoles describe induced fields only.

Nodal interpolation alone would mismatch the analytic multipole trace between
outer nodes. Each outer patch therefore includes a transfinite correction:
normal linear blend times [analytic outer trace - polynomial trace]. It vanishes
on the inner plasma and patch interfaces and exactly restores the analytic trace
on the outer matching surface. Both derivatives enter the SAME energy Gram
matrix. No unpaired boundary field substitution or hidden current sheet is
introduced at the outer matching surface. Interior normal B is continuous by
construction; tangential B matching remains a convergent variational condition.

The assembled energy is (pi/mu0) integral [(psi_r²+psi_z²)/r] dr dz plus the
analytic exterior. Vacuum internal nodes and multipoles minimize this energy
at fixed plasma boundary flux; the Schur complement is the exterior Hessian.
Diagonal scaling improves the positive solve without changing the equations.
Only roundoff antisymmetry is symmetrized. The full matrix and extension response
are retained to evaluate physical exterior fields, including the regular axis.
Finite patch resolution and multipole order are explicit approximations.

This is the required finite-domain magnetic component of frc-vacuum.2. The full
finite current/field maps, source work and moving plasma integration are still
pending; this component alone is not a finite-coil plasma simulation. Independent
dipole fields and exterior energy provide mathematical physical comparisons.
