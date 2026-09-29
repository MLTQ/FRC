# Moving-material verification fixtures and records

`MovingCase` declares all numerical controls and builds independent unit-scale
gas-expansion, initialized spherical-FRC and prescribed Gaussian-bias fixtures.
Density is 1 kg/m³ and gamma is 5/3; these are mathematical benchmarks, without
a material species or transport interpretation. Exterior magnetic spacing is
held constant when the matching sphere moves. Source amplitudes are tesla.

`errors` compares expansion against lambda=sqrt(1+2t²) and measures equilibrium
acceleration, actual surface displacement and ideal nodal-flux drift. The
Gaussian case has an initial equilibrium corresponding to its initial bias;
its subsequent motion has no claimed analytic solution. `snapshot` retains
the complete accepted material state, constants and vacuum coefficients in
lossless arrays, so records can be reconstructed without hidden initialization.

`provenance` hashes all production Python sources and three measurement/fixture
modules along with the running environment. Report/plot code is separately
hashed when assembling the completed study. Old fixed-boundary studies retain
their historical fingerprints and are not silently relabeled as current.
