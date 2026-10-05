# Explicit magnetic wall space

Wraps an unchanged compatible geometry and restricts Btheta(R)=0 with its exact
radial histopolant trace. The magnetic-mass-orthogonal projection defines the
new variational curl Fw=P F. Original field sampling, magnetic Hessian, electric
mass, poloidal flux and gradient are unchanged. P F grad=0 and poloidal div B
remain compatible. The separate original geometry is retained for controls.

`pack` accepts only admissible toroidal states, allowing normalized floating
trace roundoff <=256 machine eps; it never changes supplied coefficients.
`project_curl` is the declared variational derivative, not a state filter. The
unprojected Faraday defect must be reported separately and converge in physical
checks. Full derivation, zero boundary work and adoption gates are in
docs/STRONG_WALL_DESIGN.md. This is a bounded-wall candidate, not vacuum.
