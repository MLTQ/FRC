# Physical bridge for the regular-coordinate fluid

Maps the actual regular fluid degrees of freedom onto compatible electromagnetic
quadrature. Scalars use the element degree-p polynomial in x=r²/2 and the real
periodic trigonometric interpolant in z. Transverse velocities use degree p−1
at the first element, with the inactive axis coefficient omitted. This is the
same field definition as the independent evolution sampler, not a cell-mean
conversion. Field radial cells must subdivide fluid elements so quadrature
does not unknowingly cross a discontinuous element boundary.

`scalar` and `transverse` are linear samples; `adjoint` is their physical
volume-weighted transpose. `velocity` samples r*S_v(q/n), r*S_v(l/n), and
S(mz)/S(n). `force_rate` maps a supplied physical force density to regular
momentum rates using the true grid kinetic-energy gradient. It pairs exact
physical force work and preserves integrated axial force and angular torque.
It is not a least-squares force fit or a global energy correction. Physical
force accuracy and quadrature convergence must be independently verified.

For the transverse components, dot(q)=S_v.T W(r*f_r)/(M*V*r_node²), with
zero inactive axis slots; the angular formula is identical. Axial dot(mz) is
n_node*S.T W(f_z/n_h)/(M*V). These formulas are density dependent where the
declared physical velocity is density dependent. `RegularDensity` supplies
the same linear scalar map and its pullback to variable-current inertia.

All arrays must be finite, point and sampled densities positive. There is no
positivity floor, viscosity, vacuum treatment or magnetic evolution here.
