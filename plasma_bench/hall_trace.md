# Zero toroidal-field trace Hall geometry

`TraceHallCurl` retains compatible psi/r theta and Er spaces but includes **all**
Ez scalar radial nodes, including R. Its physical Faraday map has the last-cell
+Ez(R)/dx term and its positive axial electric mass includes that basis.
The adjoint magnetic current therefore uses the natural prescribed Btheta(R)=0
trace. Psi is fixed at the axis/wall, hence Br(R)=0 and Etheta(R)=0. Ez(R) is
free: this is a mixed magnetic boundary, not a conducting Ez=0 wall or vacuum.

The complete Q1 scalar gradient belongs to this electric space. C grad=0,
including at the wall, and its weak current pairing matches `DualCurrent` with
zero toroidal wall trace. Normal current and magnetic boundary power vanish.
Free Ez must not be combined with an extrapolated nonzero toroidal wall field
in the energy flux. `compatible_poynting` honors the explicit `wall_btheta`.

`scalar_samples` interpolates a full dual scalar with positive Q1 weights onto
the common Hall quadrature. Warm carrier density uses this explicit placement,
not a different reconstruction inferred from primal density. The legacy
conducting Hall geometry and its recorded baseline remain separate.
