# Primal-to-dual conservative mapping tests

Known uniform and Q1 volume moments check overlap conservation and its energy
adjoint. Independent affine contravariant physical fluxes have exact dual-face
integrals and constant divergence, including half axis/wall supports. Random
nonzero boundary fluxes check **local** commuting divergence at every cell and
the global boundary balance, with periodic seam duplication rejected if wrong.
Constant specific entropy follows compressing density under both reconstruction
orders; this checks an exact invariant needed for species pressure work.
