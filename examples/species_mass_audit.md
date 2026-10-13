# Physical mass and representability of the entropy source

## Purpose
Read the fully replayed weak entropy audit and split its strong-source discrepancy
into the actual interface/transfer remainder, diagonal-versus-physical mass map,
and the unresolved physical polynomial/Fourier projection remainder. No solver
or accepted state is changed. The design is SPECIES_MASS_DIAGNOSTIC_DESIGN.md.

## Components
- `mass_reference`: independent analytic normalized mass matrices, using the
  highest local Legendre mode and real Fourier Nyquist cosine. Lower modes
  retain unit mass ratios; the exceptional ratios are 5/11 and 1/2.
- `run`: independently integrates q16/20 tensor masses, solves for diagnostic
  L2 coefficients, reconstructs the complete field difference identity, and
  records pairwise products as well as norms and source modal fractions.

## Contracts
The original K and every component vector come from the prior exact replay.
The physical target is the actual independently reconstructed current acting
on S(C/n), with the original species multiplier and base density. Every species
rate uses the same SI units; no unlike-component norm is formed. Projection
error remains separate from the diagonal metric's effect. Global conservation
residuals retain absolute values and corresponding absolute-integral scales.

Exact-mass diagnostic coefficients are not passed to any production solver.
This audit cannot identify the cause of coupled growth by itself. Complete
fresh replay, q16/20 differences and analytic mass checks are required.
Run under `guarded_species_mass`; no unit tests or physical repairs.
