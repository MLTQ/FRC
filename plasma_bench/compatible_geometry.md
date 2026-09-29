# Higher-order compatible cylindrical geometry

Uses x=r²/2 and common physical r,z Gauss quadrature. N denotes nodal
interpolation, E interval-average edge histopolation and D incidence. Fields are
Br=−Nx Ez Dz psi/r, Bz=Ex Dx psi Nz, Btheta=r Ex g Ez.
Electric fields are Er=r Ex er Nz, Ez=Nx ez Ez, and Etheta=Nx(r_i et)Nz/r.
The derivative identity dN=E D preserves the existing Faraday incidence and
C grad=0. Full periodic psi nodes and cell-average g remain the coefficients.

Physical component masses and magnetic Hessian use separable 1-D Gram matrices,
avoiding large common-point sparse maps. Tensor evaluation still uses the same
physical quadrature, enabling independent energy and work checks. The actual
quadrature order is explicit; exactness or convergence must be measured.

Fixed zero axis psi and fixed wall psi, zero wall Etheta, natural Btheta(R)=0
and free Ez(R) match the bounded mixed trace. Axis psi is constrained before
sampling; its singular N_axis/r basis is excluded from magnetic energy. This is
not vacuum or arbitrary external-boundary geometry. Nodal scalar interpolation
is not positive. No thermal/density representation or nonlinear coupling is
silently switched to this geometry; it is a separate verified component.
