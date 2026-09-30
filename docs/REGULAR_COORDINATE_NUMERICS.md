# Regular-coordinate joint Euler transport

The mixed-flow study exposed a physical representation mismatch in the previous
Gauss-r prototype. Its conserved nodal energy can badly underestimate the
continuous energy of its declared regular r² interpolant. A new formulation
addresses the representation and axis metric directly; the earlier source and
all comparison data remain unchanged. This is still a separate fluid component,
not a replacement of the coupled formation solver.

## Actual state and original equations

Let x=r²/2, u=vr/r and Omega=vtheta/r. The stored point variables are
U=(n,q=n u,l=n Omega,mz=n vz,Ce=n se,Ci=n si). Actual physical angular
momentum is L=2x l, **not l**. Scalars and regular velocities are finite at x=0.
The original axisymmetric Euler equations become

* n_t+∂x(2x n u)+∂z(n vz)=0;
* q_t+∂x(2x q u)+∂z(q vz)=n(Omega²−u²)−p_x/M;
* l_t+∂x(2x l u)+∂z(l vz)=−2n u Omega;
* mz_t+∂x(2x mz u)+∂z(mz vz)=−p_z/M;
* Ca_t+∂x(2x Ca u)+∂z(Ca vz)=0.

Here M=mi+me and p=pe+pi. The physical two-species entropy/EOS convention is
unchanged: ha=e exp(sa)(n/n0)^(gamma_a−1)/(gamma_a−1), pa=(gamma_a−1)n ha.
The manufactured gamma choices test unequal species laws; they do not supply
physical ionization/collision closure for the eventual formation experiment.

Physical energy is H=M x(q²+l²)/n+M mz²/(2n)+n sum ha. Its true gradients are

wn=−M x(u²+Omega²)−M vz²/2+sum(gamma_a−sa)ha;
wq=2M x u, wl=2M x Omega, wmz=M vz, wCa=ha.

## Flux and geometric work

Use positive logarithmic means and arithmetic bars. For a radial pair define
speed a_bar=x_i u_i+x_j u_j=bar(2x u), Gamma=n_log*a_bar. The fluxes are

Fn=Gamma; Fq=Gamma*bar(u)+bar(p)/M;
Fl=Gamma*bar(Omega); Fmz=Gamma*bar(vz);
FCa=a_bar*[n_log*bar(sa)+(gamma_a−1)(bar(n)−n_log)
                 +n_log*(bar(ha)/h_log−1)].

The scalar energy potential is psi=2x p u. Direct original-EOS algebra gives

[w] dot F = [psi]+M Gamma [x]*(u_i u_j+Omega_i Omega_j).

For Q=W D, the integrated geometric terms are

Sq_i=sum_j Qij[x]n_log*(Omega_i²−u_i u_j);
Sl_i=−2u_i sum_j Qij[x]n_log*bar(Omega).

Their true kinetic work cancels the extra pair jump. The weighted angular
identity also holds: flux differencing of 2x*l has an extra term canceled by
2x*Sl. Thus n,L,mz,Ce,Ci and actual H conserve in the declared quadrature.
Neither q nor the unweighted l sum is a global invariant. The pressure term
is evaluated as sum Qij(p_j−p_i)/M to avoid constant-pressure cancellation.

Equal-x axial pairs use the same thermal bracket with speed bar(vz), fluxes
Gamma*(u,Omega,vz), and pressure only in axial momentum. Periodic skew axial
flux differencing remains order2/4/6. No global energy or force correction is
introduced. RK4 finite-step drift is measured separately.

## Dependent axis velocities and physical representation

Each radial element uses p+1 Gauss-Lobatto points in x and positive quadrature
weights, including its endpoints. The first element has no independent q/l
velocity at x=0: u and Omega there are traces of degree-p−1 polynomials through
the p positive nodes. Two packed axis slots stay zero and are explicitly
inactive. Other scalars, including density and entropy, remain active there.
This is a DOF definition, not a post-step energy or temperature repair.

The reason matters: an independent velocity at x=0 would have zero nodal
kinetic weight while affecting a continuous interpolant. Removing that degree
of freedom gives a positive physical velocity metric. For constant density,
the first-element integral of x*u² is exact, since its polynomial degree is
2p−1, precisely the Lobatto integration degree. The analytic f(r)=r work
adjoint then gives the correct regular acceleration u_t=constant.

Continuous density,mz,Ce,Ci are degree-p polynomials in x. u/Omega are
degree-p−1 in the first element and degree-p elsewhere, through active q/n and
l/n values. Physical momentum is n_h*r*u_h and L=n_h*r²*Omega_h. Axial
velocity is mz_h/n_h. All active point values are retained; no filtering or
fit is used. The field is regular at exactly r=0. Variable-density energy and
angular integrals are not generally integrated exactly; their discrepancy
from grid quantities is saved independently and must converge.

Shared centered interface fluxes and pressure-only boundary momentum fluxes
give zero particle/energy/angular flux at the exact axis and outer wall.
There is no positive radius floor. The domain remains fixed-radius, periodic
in z, smooth and positively filled. No vacuum, diffusion, limiter or shock
entropy inequality is supplied. Sampled density/pressure must remain finite
and positive, and rejected stages preserve the accepted state.

## Coupling gate

Current/magnetic coupling must use this actual density and velocity definition,
including the density pullback of current kinetic energy and actual angular
torque. Common physical force tests must precede adoption; an energy identity
alone cannot justify a mismatched force adjoint. Full drift/Maxwell/species
stress and enthalpy/Poynting work, physical transport and ionization, exact
vacuum outside a heated moving plasma, real ends/coils and converged formation
remain open. No 3-D stability or particle-confinement claim follows.
