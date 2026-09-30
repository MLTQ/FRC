# Cylindrical joint physical transport

This separate smooth, positive-density axisymmetric fluid prototype extends
the axial joint flux. It has no magnetic field, current, free boundary or
formation experiment. It has not replaced the coupled plasma solver.

## Original equations and actual state

Let M=mi+me, p=pe+pi, D=∂r(r vr)/r+∂z vz and material derivative
Dt=∂t+vr ∂r+vz ∂z. The verification target is

* Dt n=−n D;
* Dt vr=vtheta²/r−∂r p/(M n);
* Dt vtheta=−vr vtheta/r;
* Dt vz=−∂z p/(M n);
* Dt pa=−gamma_a pa D, separately for each species.

The actual point state is U=(n,mr=n vr,L=n r vtheta,mz=n vz,Ce=n se,Ci=n si).
There is no fitted or separately repaired angular momentum. Entropy convention
is sa=log[pa/(n e (n/n0)^(gamma_a−1))], n0=1e20/m³, with temperature in eV.
Physical h_a=e exp(sa)(n/n0)^(gamma_a−1)/(gamma_a−1), and

H=M[mr²+(L/r)²+mz²]/(2n)+n sum_a h_a.

Energy variables are wn=−M|v|²/2+sum_a(gamma_a−sa)h_a,
wmr=M vr, wL=M vtheta/r=M Omega, wmz=M vz, wCa=h_a.
The constant species masses and ideal-gas gammas are fixture inputs. Unequal
gammas stress the two-species implementation; they do not establish a physical
argon ionization/transport closure. The pressure laws contain no dissipation.

## Radial geometric work

Use the symmetric axial two-entropy flux of TWO_ENTROPY_NUMERICS for physical
momenta, but set the radial angular flux to

F_L=bar(r) Gamma bar(vtheta), Gamma=n_log bar(vr).

For a pair at different radii, with jump [a]=a_j−a_i and psi=p vr,

[w] dot F = [psi]−M Gamma [r] bar(vtheta) bar(Omega).

This extra kinetic work must be paired with the centrifugal force; it is not
canceled by changing the conservative angular flux or by repairing energy.
With radial area aij=bar(r), the integrated centrifugal source at node i is

S_mr,i=sum_j QN_ij aij n_log [r] bar(vtheta) bar(Omega).

Its work M vr_i S_mr,i cancels the pair's extra angular kinetic work after
summation. The pressure flux and geometric pressure source combine as

−sum_j QN_ij aij (p_j−p_i)/M.

This is algebraically the conservative r p flux plus the cylindrical p source,
using QN 1=0. Evaluating a pressure difference avoids cancellation of a large
constant pressure against its geometric source. It is not a fitted force.

## Gauss volume, boundary traces and conservation

Each radial element has p+1 Gauss nodes in r, W quadrature, polynomial derivative
D and endpoint interpolation R. Q=W D, B=diag(−1,+1), and

QN=[[Q−RᵀBR/2,RᵀB/2],[−BR/2,B/2]], lift=[I,Rᵀ].

QN 1=0 and QN+QNᵀ=diag(0volume,B). The volume/boundary construction follows
[Chan's quadrature-based SBP framework](https://arxiv.org/abs/1708.01243).
Our two-species physical-energy variables and cylindrical source pairing are
derived here; the citation is not a claim that the paper supplies an FRC model.

Boundary traces interpolate **w**, then invert the same physical EOS at the
actual boundary radius. The inversion gives vtheta=r Omega and L=n r² Omega
even at exactly r=0. No volume node sits at zero and there is no radius floor.
Normal velocity is imposed weakly through zero mass flux at axis/outer wall;
the boundary pressure-only momentum flux gives zero physical energy flux.
Adjacent elements share the same symmetric interface flux. Axial derivatives
use periodic order2/4/6 skew flux differencing at equal radius.

Conserved integrated grid quantities are n, L, mz, Ce, Ci and **actual H(U)**
weighted by 2 pi r dr dz. Radial momentum is not a cylindrical invariant.
Energy conservation uses the interpolated true H gradient and shared interface
energy flux G=bar(w) dot F−bar(psi). The independent chain-rule check balances
each radial element's energy, not only a global ledger. Finite-step RK4 energy
drift is separately measured; there is no evolved energy account or reset.

## Regular continuous physical interpretation and limits

For independent integration, reconstruct n,mz,Ce,Ci,mr/r,L/r² as piecewise
polynomials in r² through the actual point values, with periodic trigonometric
interpolation in z. This gives even scalars, odd radial/azimuthal velocities
and finite physical kinetic energy at the axis. A naive arbitrary polynomial
L(r) can have nonzero L(0) and a divergent kinetic integral; it is not used.

This regular physical reconstruction differs from the polynomial in r used by
SBP differentiation. Nodal energy/angular conservation does not imply exact
conservation of its continuous energy/angular integral. Both differences and
higher-quadrature convergence are retained. Regular interpolation also does
not prove that the evolving numerical solution is accurate at the axis.
Original PDE first-element/max norms and radial mode controls remain required.

The driver adds no viscosity, heat, floor, filter, post-step repair or atmosphere.
All stages and boundary thermal traces must be positive/finite; rejection is
atomic. This nondissipative method is not shock capturing and does not guarantee
positivity or nonlinear stability. A common representation/work bridge to the
existing current and magnetic fields remains unproved. Exact vacuum outside a
heated moving plasma and physical transport/ionization still require separate
work. Neither this axisymmetric fixture nor closed flux proves 3-D stability
or particle confinement.
