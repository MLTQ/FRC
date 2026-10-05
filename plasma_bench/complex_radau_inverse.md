# Real two-stage Radau inverse

For A=[[5/12,-1/12],[3/4,1/4]], the conjugate eigenvalues are
1/3 +/- i sqrt(2)/6 with eigenvectors ((1 +/- 2 i sqrt(2))/9,1).
Transform real residuals to one conjugate component, apply the existing
complex physical block inverse at dt times that eigenvalue, and reconstruct
the two real stage corrections. Verify by substitution into the real stage
equation. This changes only the preconditioner, never a physical RHS term.
