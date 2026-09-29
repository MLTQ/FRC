# Independent compatible-field norms

Evaluates the declared nodal/interval-average polynomial representations using
root-product Lagrange polynomials and physical cylindrical quadrature, without
production basis solves, sampling, energy matrices or current reconstruction.
Reads saved full B/J/Bdot/Jdot/full-E cochains and compares to independent
analytic parent-wave fields/derivatives. Reports physical and first-annulus L2,
dimensional RMS, sampled maxima, magnetic energy and physical reduced-mass
current energy. A prescribed bias is subtracted only in the perturbation norm;
full magnetic energy and saved inputs retain it.

Affine records explicitly supply `axial_bias_t` beside increment cochains.
Physical energy adds that bias before squaring, preserving local cross energy.
Perturbation norms use the independently stored increment directly, avoiding an
unnecessary add/subtract loss of digits. Older rounded-total records retain their
original interpretation; their actual loss of precision is not repaired by the
audit. Both representations and all denominators remain explicit in metadata.

`direct_basis` is a separate polynomial algorithm for the explicitly defined
local representation. It does not imply that sharing a mathematical field
definition is an independent test of its stability. Different quadrature orders,
polynomial exactness, original mixed equations and physical wave convergence
provide complementary verification. No state/trajectory is changed.
