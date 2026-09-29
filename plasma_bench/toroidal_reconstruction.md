# Toroidal reconstruction and cylindrical stress

Fits even variables and regular Omega=vtheta/r, g=Btheta/r in x=r²/2. For
second order, vr/Br use the established odd volume moments; first order keeps
piecewise-constant poloidal states. Toroidal total pressure includes exact
<r²> g²/2mu0. Radial traces use actual radius; axial quadratic toroidal moments
use sqrt(<r²>). The scheme is exact for the rigid-rotation/linear-pinch fixture
and second-order consistent for smooth variable subcell profiles.

Shared MC/positivity helpers limit reconstructed density and thermal traces;
no conserved average changes. Every face, corner and source-quadrature trace
participates in the pressure blend. The source integrates the same reconstructed
pt+rho r² Omega²-r² g²/mu0 with five-point Gauss (degree eight for the rotation
term). This includes centrifugal and magnetic hoop stress with correct signs.
Statistics are exposed for accepted-stage accounting. Positive traces do not
prove positivity of the full update. See TOROIDAL_MHD_NUMERICS.
