# Regular nodal bulk velocity

Q1 interpolation uses x=r²/2 and periodic z on the complete dual scalar nodes.
scalar returns value, partial_x and partial_z at the common physical Hall
quadrature. velocity accepts explicitly regular components (vr/r,vtheta/r,vz),
so finite axis values are supplied rather than guessed by dividing zero velocity
by zero radius. It returns physical cylindrical v and curl(v), including the
azimuthal-basis derivative 2 Omega+r² partial_x Omega. No finite-difference
current derivative, axis floor or hidden boundary extension is used.

This defines a representation, not a primal-to-dual momentum map. A future
moving fluid coupling must supply a conservative, independently verified map.
Scalar density sampling agrees with TraceHallCurl; no extra density is evolved.
