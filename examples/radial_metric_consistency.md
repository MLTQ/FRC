# Radial polynomial and metric consistency probe

Compares the actual cell-momentum reconstruction R with the continuity-lifted
physical momentum L=Q B. Q maps independent interior radial particle faces to
physical momentum; B maps actual cell momentum means to those faces. For an
analytic regular force r/R(1−r²/R²)², which also vanishes at the wall, compares
the old mixed work projection (R.T W R)^-1 L.T W f with the consistent face
projection (Q.T W Q)^-1 Q.T W f. The target is an independent polynomial, not
a force taken from any numerical solver or equilibrium cancellation.

Records direct representations, force errors including the first annulus,
rank/nullspace evidence and both physical Gram matrices across degree3/5 and
nr8/16/32/64. It neither changes the frozen solver nor adopts a new momentum
state. A face-state fluid would require its own true energy, convection,
conservation, parent-rate and trajectory verification before use.
