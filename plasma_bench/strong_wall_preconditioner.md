# Wall-restricted approximate inverse

Constructs the frozen original uniform-background block operator, then restricts
it to the admissible wall coordinates. The magnetic left restriction uses
K=H+A.T I(n) A with the *uniform preconditioner* density and physical masses.
The same helper can form K at a supplied nonuniform diagnostic state. Fluid
test/trial coordinates omit normal wall momentum. The actual nonlinear RHS
uses the complete StrongWallPlasma equations at the evolving physical state.

`magnetic_metric` builds the full scaled magnetic/current Hessian;
`magnetic_left` returns (Q.T K Q)^-1 Q.T K. Both are dense bounded diagnostic
operations. StrongWallPreconditioner inherits the original block factor/solve;
ComplexPlasmaInverse consumes the same four blocks. It supplies proposals only,
never an accepted derivative, changed mass, damping or energy correction.
