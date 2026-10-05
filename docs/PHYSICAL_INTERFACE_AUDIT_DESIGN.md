# Original physical fluxes and one-sided interface audit

This is an independent audit design for the existing bounded, periodic,
collisionless physical-mass model. It changes no evolution equation or state.
The remaining interface gate cannot be replaced by a small global energy drift.

Let M=mi+me, alpha=me/M, beta=mi/M, mu=mi*me/M, w=J/(e*n), with
vi=v+alpha*w and ve=v-beta*w. The original species parent energy flux is

    F = sum_s vs * [ms*n*|vs|^2/2 + gamma_s*Us] + E cross B / mu0.

The density is shared and the continuum requires div J=0. The physical energy
is M*n*|v|^2/2 + sum Us + mu*|J|^2/(2*e^2*n) + |B|^2/(2*mu0).
For conservative C_s=n*s_s, its thermal rate is
h_s*[C_s,t + (gamma_s-s_s)*n_t], using the actual complete numerical RHS.
Current and magnetic rates use actual J_t and B_t; constant harmonic bias is
included in B dot B_t and E cross B. No unresolved term is subtracted or repaired.

Each field cell is also aligned with fluid radial-element boundaries. Integrate
actual energy/particle/entropy time derivatives in its volume, and integrate
original physical fluxes on each face using the trace *inside that cell*.
Particle flux is n*v, charge-current flux is J, entropy flux is C_s*vs.
A normal-current jump is measured, not set to zero by an assumed quasineutral
identity. Spectral fluid dependence is periodic; radial fluid states may have
different one-sided element traces. Field reconstruction has continuous
appropriate magnetic normal/electric tangential traces, but projected current
normal traces need independent inspection.

At each internal face retain both traces, their integrated jump and the
quadrature integral of the absolute jump. Do not use small net cancellation
as proof of a small local defect. The periodic axial seam is an internal face.
The axis has exactly zero cylindrical face area; the outer wall uses the
interior trace. The sum of all cell outward fluxes equals exterior flux plus
left-minus-right internal traces. This distinguishes classical piecewise
volume divergence from its distributional jump contribution.

Evaluate traces at the exact face with an explicitly selected polynomial cell,
not by displacing points with an epsilon or silently averaging two states.
Independent root-product reconstruction and geometric Gauss face quadrature
must be tested against known smooth and discontinuous manufactured functions.
Then measure complete full-RHS local budgets and their mesh/quadrature behavior,
including separate total-energy, particle, current and species entropy results.
These measurements do not themselves establish a local conservative repair,
heated moving vacuum, laboratory wall/sheath physics or formation readiness.
