# Full-vector moving finite-mass current and plasma

Purpose: evolve poloidal and toroidal magnetic fields, all carrier-current
components, bulk rotation and meridional motion with both species and heat.
This extends the physical material equations rather than freezing poloidal flux
or adding a diagnostic current. A complete coarse trajectory and its independent
force/virial comparison are recorded in docs/MATERIAL_VECTOR_COUPLING.md; the
higher-order full trajectory is pending. Coarse local fields are not converged.

The state and unsplit RK4 remain MaterialCurrentState/MaterialCarrierCurrent.
MaterialVectorGeometry supplies full Piola maps, including angular winding and
metric derivatives. Ampere is e*Mn*w=-D.T*b and Faraday H*bdot=D*E*. The full
force in Ohm is delta*w cross (B+mu/e curl w)-mu/e*w cross curl v, plus pressure
and eta*J. All physical mass factors and current kinetic energy remain.

Magnetic virtual work is the full Maxwell stress. Current virtual work is
mu*n*w*w:grad(delta x). Both include the toroidal curvature force and angular
torque. The latter is -integral r*Btheta*(Bpol.grad N)/mu0 plus
integral mu*n*r*wtheta*(wpol.grad N). No rotation or torque is prescribed.
Their mechanical work uses actual nodal velocity and angular velocity and
cancels the magnetic/current metric work. Axis radial motion alone is constrained.

MaterialDualEntropy provides conservative upstream species flux and its matched
pressure work. Its first-order numerical diffusion remains explicit. Toroidal
carrier flow has no entropy derivative for axisymmetric scalars. Prescribed
constant resistance supplies actual J² heat with the mass-weighted species share;
optional conditional Maxwellian exchange uses both evolving temperatures.
No heat suppression, atmosphere, density/entropy floor or energy reset is used.

The physical mixed-row residual tolerance remains 1e-10. Initial projection
minimizes physical vector B and carrier-drift errors against consistent analytic
targets. It does not fit evolved states. The physical laboratory E restores
thermal/self-drift potentials and subtracts v cross B from material E*.

The present vector space has B=0 on the moving side/ends and exactly zero field,
matter and conductivity outside. Penetrating vacuum fields and finite coils are
still absent. Initialized closed psi is not resistive formation. Cold-edge
transport, pointwise current traces, spatial convergence, 3-D stability and
particle confinement remain unproven even if energy and torque balances agree.
