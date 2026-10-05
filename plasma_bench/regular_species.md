# Compatible regular species pressure and entropy

## Purpose
An isolated bounded instantaneous operator on the existing actual regular state.
It pairs physical electron/ion pressure work with their entropy-density rates.
It is not itself a complete moving integrator or a vacuum/formation model.
New RegularPlasma composes it with fluid/current/magnetic evolution and verifies
bounded complete fast-wave periods; see REGULAR_PLASMA_VERIFICATION for scope.

## Components
`SpeciesPotential` samples n*f and n at field nodes, averages both conservative
traces at radial interfaces, and returns T f=N(n f)/N n. Its adjoint is the
exact Euclidean transpose. Positive point, quadrature and common-trace densities
are required; no floor is applied. Axial interpolation is the existing periodic
Fourier interpolant, hence the graph is not compact in z.

`RegularSpecies` uses the physical L2 `SpeciesProjection` by default. Explicit
`transfer='point'` retains the handoff trace map as a failed accuracy control.
`RegularSpecies.transport` applies K=(A-A^T)/2, with
A=S^T W J_h dot(P G T), using tensor actions without a dense matrix.
`skew=False` preserves the handoff proposal as an explicit diagnostic control.
`pair` returns Cdot_s=-a_s K s/(e V), zero other fluid components, and
F_s=a_s/e [M_E G T(gamma*h-s*h/2)
-P^T W((S h)P G T s-(S s)P G T h)/2]. Here a_e=-beta, a_i=alpha.
The changed potential is required by the skew pressure adjoint.
The code stores an exact affine split with ha=h[0], sa=s[0], dh=h-ha, ds=s-sa.
The reduced term uses dh*grad(ds)-ds*grad(dh), and the potential additionally
contains -ha*s/2+sa*h/2. Constant h or s then stays in the curl kernel without
subtracting nearly equal full fields. These are algebraic anchors, not a fitted
equilibrium or modified temperature. `reduced_forcing` drives induction;
`potential_electric` is restored only in the full electric diagnostic/equation.

## Contracts
The actual EOS derivative h, common point density and physical masses are kept.
Inputs are evaluated in float64 without mutating the caller; integer or float32
arrays cannot truncate the returned physical rates to their original dtype.
For an Ampere current, A*1=0 and 1^T*A=0 by the compatible curl/gradient identity
up to roundoff. Thus K*1=0, sum V*Cdot=0, and frozen-density sum V*n*s^2/2
is constant. J dot F equals sum V*h*Cdot; no global correction is subtracted.
Pairwise graph charge flux 2 K_ij is antisymmetric and has zero divergence;
its entropy flux is K_ij(s_i+s_j). This graph interpretation includes shared
radial traces but does not establish continuum local flux accuracy.

## Limits
Zero numerical exterior current flux and periodic z match CompatibleGeometry.
The physical reconstructed current need not have zero strong divergence or
normal jumps. No numerical diffusion, heating, positivity limiter or temporal
method is included. Frozen-density quadratic stability is not nonlinear
thermodynamic, shock, coupled-plasma or formation stability. Local rate accuracy,
continuous energy and interface consistency must be measured independently.
