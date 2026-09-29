# Smooth cubic reconstruction of annular moments

radial_matrix fits four actual annular cell averages to r^s P3(x), x=r²/2,
with s=0 for even scalars, s=1 for radial momentum and s=2 for angular momentum.
Five-point physical quadrature integrates every basis moment exactly. A local
shift and scale keep the moment system conditioned near the axis and on fine
meshes. The interpolant uses interior stencils only, without invented exterior
matter. axial_samples reconstructs periodic volume averages with a cubic.

normal_momentum_flux is a candidate smooth continuity flux for the explicitly
impermeable cylinder. It directly interpolates momentum, since continuity flux
is momentum, and retains exact zero radial boundary flux. It is currently a
diagnostic only: FiniteMassMHD does not use it. Changing continuity alone would
leave bulk enthalpy/kinetic/momentum work inconsistent, so full coupled transport
must be derived and tested before adoption.

There is no slope limiter, numerical diffusion or positivity guarantee in this
interpolator. Arbitrary input cell averages are not repaired; sharp profiles
can overshoot. Polynomial exactness, smooth convergence and an independent
nearly solenoidal wave diagnose its intended role without claiming shock or
formation robustness. All inputs remain unchanged.
