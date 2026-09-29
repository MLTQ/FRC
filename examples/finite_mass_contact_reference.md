# Independent finite-mass translating contact

Exact unforced parent species share constant axial velocity with arbitrary
translating positive n and complementary pressures pe+pi=constant. The explicit
mass density is (mi+me)n; gammas can differ. Uniform axial B has J=0 and full
Ez=-pe_z/(ne)=pi_z/(ne), which is curl free and balances both species momentum.
Pressure equations are exact advection since div v=0, and full total energy
includes the physical total mass. No constraint force or manufactured heat is used.

fields(rate=True) returns analytic time derivatives of density, both pressures,
momentum and total energy. The implementation imports no production solver or
geometry. Finite amplitudes test nonlinear entropy reconstruction and changing
density; a small-amplitude wave alone would not establish this property.
