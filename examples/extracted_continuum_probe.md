# Matched continuum action on the extracted response

Loads hash-checked normalized conservative direction and actual full-RHS tangent
from the extracted-response probe. Independently reconstructs both in physical
space, then projects each onto the previously verified Fourier/Bessel parent
operator. No production matrix or accepted state is changed. Compares the
projected actual rate with the independently derived continuum action on the
same projected direction, component by component with physical radial norms.

Reports captured direction and rate fractions, radial-mode power distribution,
projected growth, actual/continuum discrepancy, and strong curl(B)/mu0 versus
weak current. Controls axial32/64, radial8/16/32 and quadrature16/20. Uncaptured
norm and strong/weak current mismatch prevent an incomplete projection from
being presented as a fully matched eigenmode. The projection has the explicit
closed-wall Bessel traces and omits the even axial Nyquist coefficient.

Saves all projected coefficients/actions/normalizations, hashes and source.
Fresh replay must agree bitwise. A discrepancy diagnoses consistency of this
specific direction; it does not prove a unique repair, formation or confinement.
