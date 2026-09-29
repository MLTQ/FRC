# Pressure-only principal symbol

Derives the local high-wavenumber symbol directly from the continuum pressure
Ohm term and electron adiabatic equation, independent of production grids.
For a unit direction k=(kr,kz), a=grad log n, c=grad pe-gamma pe a,

    delta(Btheta, pe)_t = i |k| [[0,A],[C,D]] delta(Btheta,pe)
    A=(ar kz-az kr)/(n e)
    C=(kr cz-kz cr)/(mu0 n e)
    D=J0 dot k/(n e).

Axis-curvature and background-gradient terms without perturbation derivatives
are lower order in this local symbol. The discriminant D²+4AC is independent
of the units chosen for the two perturbation amplitudes. A negative value
gives exponential growth |k| sqrt(-discriminant)/2 in this frozen subproblem.
The actual background current advection is retained; it is not set to zero.

Scans an explicitly saved physical position/direction grid of the original
ElectronReference and records one counterexample, its matrix/eigenvalues and
source hashes. Uniform density supplies a zero-A control. This is a local
necessary-hyperbolicity diagnostic, not a nonlinear well-posedness proof for
the full Hall system or a fit to the discrete failure. Do not attribute unrelated
incompressible Hall ill-posedness theorems to this thermal splitting mechanism.
