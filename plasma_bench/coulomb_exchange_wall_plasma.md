# Electron-ion temperature exchange in moving plasma

Adds thermal equilibration to the full moving finite-mass wall model. It
requires Z=1 (the existing equal electron/ion density state), both gamma=5/3,
weakly coupled Maxwellian species and small relative drift. Scalar resistivity
and its heat fraction remain explicit, separately prescribed inputs; this does
not silently install a full Coulomb conductivity tensor.

For temperatures in eV and density in m^-3, define C=e^4/(4 pi epsilon0)^2:

    nu_E = 8 sqrt(2 pi) n C lnLambda /
           [3 me mi (e Te/me + e Ti/mi)^(3/2)]
    Qe = (3/2) n e nu_E (Ti-Te), Qi = -Qe.

This is the SI Maxwellian Landau equilibration coefficient. The
[NRL Plasma Formulary 2023, pp.33–34](https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf)
gives the rounded cgs/eV coefficient1.8e-19. Its low/high-Te Coulomb-log branches
are evaluated on each actual state. Unsupported slow-ion ordering, nonpositive
log or nonfinite rate rejects explicitly. Small positive logs are reported for
review, not clipped. This is a conditional closure, not a validity certificate;
large drift, strong coupling, non-Maxwellian distributions, neutral collisions,
ionization, radiation and thermal conduction require additional physics.

The source is collocated on actual fluid nodes. With h_s=dU_s/dC_s, add
Cdot_s=Q_s/h_s, leaving mass and momentum unchanged. Therefore each node has
zero net thermal work, and its physical entropy production is
Qe/e*(1/Te-1/Ti)>=0. No global energy reset, floor or heat clipping is used.
The collocation is a spatial approximation; continuously reconstructed source,
energy and entropy need convergence measurements for nonuniform cases. The
uniform relaxation example exercises complete moving/field evolution, while
providing an independent scalar temperature reference.

The original Gauss stages see the full augmented RHS. Existing fixed psi,
zero wall Btheta/normal velocity, periodic axial domain and regular axis remain.
There is no exact-vacuum interface or formation model in this filled cylinder.
No numerical damping is added by this class. It does not combine with the
separate hyperdiffusive candidate by implicit or hidden inheritance.

Two complete400us uniform relaxation runs at40/80 steps show the expected
time refinement against an independent scalar integral. See
docs/COULOMB_EXCHANGE_RESULTS.md. Nonuniform advection/source coupling is the
next physical simulation (frc-transport.1), not established by the uniform case.
