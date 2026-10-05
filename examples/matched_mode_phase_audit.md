# Representative fast-frequency phase diagnostic

Uses the independently tested uniform-rest modal physical operator on each
actual mesh and fixed mean physical scales, retaining all Fourier/radial modes.
Saves every complex eigenvalue and the uniform physical state. This is a
representative linear operator, not the actual nonuniform Jacobian, an upper
spectral bound, a projected perturbation amplitude or a physical instability.

For two-stage Gauss, R(i*x)=(1+i*x/2-x²/12)/(1-i*x/2-x²/12). Its continuous
phase from x=0 is 2*atan2(x/2,1-x²/12), tending to2*pi for large x. Reports omega*dt,
phase-frequency errors and exact uniform modal endpoint error at the actual
30/60/120-step controls. Unit-modulus stability does not imply resolved phase.
No mode is damped, removed or substituted into the evolving nonlinear model.
The results motivate time/mesh controls, not a declaration of their success.

Complete source/runtime and spectra freeze and replay bitwise. Roundoff real
parts are reported relative to the frequency scale rather than called growth.
