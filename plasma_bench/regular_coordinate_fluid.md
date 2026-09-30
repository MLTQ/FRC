# Experimental regular-coordinate Euler transport

This is a new prototype, not an adoption or mutation of the frozen cylindrical
solver. x=r²/2, u=vr/r, Omega=vtheta/r, and U=(n,q=n u,l=n Omega,mz,Ce,Ci).
Actual angular momentum is L=2x l. Physical kinetic energy is
M x(q²+l²)/n+M mz²/(2n). The original regular equations are
q_t+div(q v)=n(Omega²−u²)−p_x/M and l_t+div(l v)=−2n u Omega.

The radial pair uses Gamma=n_log*(x_i u_i+x_j u_j), q flux Gamma*bar(u)+
bar(p)/M, l flux Gamma*bar(Omega), and the same thermal entropy bracket with
transport speed bar(2x u). For true H gradients the extra kinetic jump is
M Gamma [x]*(u_i u_j+Omega_i Omega_j). Sources
S_q,i=sum Q_ij[x]n_log*(Omega_i²−u_i u_j) and
S_l,i=−2u_i sum Q_ij[x]n_log*bar(Omega) pair this work. The l source also
gives conservative **weighted actual angular momentum** 2x l. This derivation
must be verified against independent original equations before adoption.

Gauss-Lobatto nodes include exactly x=0. There are no independent radial or
angular velocity DOFs there: their trace is the degree-p−1 polynomial through
the first element's positive velocity nodes. Two stored axis slots stay zero
and are explicitly inactive. `fields` reconstructs physical traces from active
DOFs; it does not repair an evolved energy or temperature. For constant density,
quadrature integrates the first-element x-weighted velocity norm exactly.
Scalars and axial velocity retain p+1 nodes, including the axis.

Shared interface fluxes and pressure-only boundary fluxes impose impermeable
axis/wall; z is periodic. The radial pressure derivative is evaluated using
pressure differences. No radius/density floor, diffusion, heat, filter or energy
reset is added. RK4 stages reject atomically. Continuous reconstruction, full
original rates, linear spectra and mixed nonlinear trajectories are still
required; no magnetic/current, free-edge or formation support is claimed.
