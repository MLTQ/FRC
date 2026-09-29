# Independent physical state replay

Recomputes quarter-mass overlap, physical cylindrical B/J and bulk velocity,
species energies and logarithmic entropy by physical r,z quadrature. It does
not use production energy matrices/sampling or production total-energy sums to
evaluate component integrals. The input current is the declared weak Ampere
representative; analytic rate/wave comparisons independently test its accuracy.

Electron density/pressure are Q1 dual fields, ion density/pressure and radial/
axial velocity use the stated primal cell representation, rotation uses r Omega.
The ion pressure comes from the declared conservative partition; this replay
checks representation/accounting, not the physical correctness of that partition.
Both species particle integrals must agree. Physical entropy remains distinct
from nodal entropy diagnostics and has no asserted monotonicity theorem.
