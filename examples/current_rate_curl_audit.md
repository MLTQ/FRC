# Physical current-rate and magnetic-rate curl comparison

For each complete accepted signed pair, independently reconstruct the centered
current rate J and strong curl of the centered magnetic rate S=curl(Bt)/mu0 at
all original quadrature16/20 points. Original J, reference and weight hashes
must reproduce. Raw coefficient/source/reference identities are checked first.

The identity J-Jref = (S-Jref)+(J-S) separates magnetic-rate derivative error
from weak/strong current reconstruction differences. It reports squared norms,
their signed cross term, full energy-decomposition defect and each component's
own continuum-normalized L2 norm. Exact radial-row integrals are retained for
all terms. The derivative is independently constructed with root polynomials,
not with the production weak curl, mass matrix or field sampler.

Strong curl is defined inside polynomial cells and excludes distributional
face jumps. Consequently J-S is not itself a discrete Ampere violation or
proof that the weak operator is wrong. Interface and boundary conditions need
separate convergence controls. No state is changed and no term removed. Every
case must have complete original strict acceptance, independent full physical
fingerprints, quadrature controls and fresh full replay. No unit tests.
