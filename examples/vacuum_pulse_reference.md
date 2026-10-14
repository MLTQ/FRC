# Independent vacuum-driven cylindrical limits

The uniform source fixes the physical exterior Bz, not the enclosed plasma flux.
For constant density/resistivity and infinitesimal amplitude,
(1-de² Delta) B_ind,t = D Delta B_ind - B_source,t, with B_ind(R)=0.
Here Delta=drr+dr/r, D=eta/mu0 and de²=mu/(mu0 n e²). The J0 eigenvalues are its
zeros qR (the old fixed-flux comparison used J1 zeros). Coefficients of a uniform
field are c_q=2/(qR J1(qR)); b_q=-c_q convolution(B_source,t)/(1+de²q²).

A polynomial lift a(r²-R²), de² a_t+D a=B_source,t/4, restores the physical
wall-current trace and accelerates the remaining series. Bz=B_source+lift+
residual J0, Jtheta=-2ar/mu0+residual q J1/mu0. The analytic time convolution is
reused from the independent reference, not a production solver or stage record.

For a nonconstant axial seed, psi=A r J1(qr) cos(kz)/q with k=2pi/L.
Vacuum matching gives q J0(qR)+k K0(kR)/K1(kR) J1(qR)=0. The exterior is the
continuum decaying r K1(kr) solution, independently of production's polynomial
axial Galerkin eigenvalues. The field/current decay rate is
D(q²+k²)/(1+de²(q²+k²)). This seed exercises nonzero vacuum field and energy.
Superposition with the pulse is a linear mathematical baseline; the simulation
retains nonlinear fluid response, so finite-amplitude differences are possible.

Bessel orthogonality/integrals: https://dlmf.nist.gov/10.22 ; modified derivatives:
https://dlmf.nist.gov/10.29 . General magnetic diffusion method:
https://web.mit.edu/6.013_book/www/chapter10/10.6.html . The finite-inertia and
plasma/vacuum matching formulas above are derived here, not claims that these
sources validate the full plasma model. No experimental comparison yet.

Completed runs, source closures, error norms and measured limitations are in
docs/VACUUM_CYLINDER_RESULTS.md. The axial seed is initialized as a mathematical
mode; any initial closed flux is not evidence of formation.
