# Compatible species pressure and entropy on the regular state

This is a bounded collisionless adiabatic **component**, on the actual state
`n,q=n*vr/r,l=n*vtheta/r,mz,Ce=n*se,Ci=n*si`. It keeps the real electron/ion
masses, common density and existing bulk continuity. There are no collisions,
ionization, neutrals, radiation, artificial atmosphere, heat, filter or energy
reset. Its prescribed-current evolution is a verification fixture, not a moving
plasma or a formation experiment. See REGULAR_SPECIES_VERIFICATION.md for data.

## Parent equations and physical thermal derivative

Let M=mi+me, alpha=me/M, beta=mi/M, a_i=alpha and a_e=-beta. Then
v_s=v+a_s*J/(e*n), C_s=n*s_s and

    C_s,t + div(C_s*v_s) = 0,
    E_pressure,s = a_s*grad(p_s)/(e*n),
    h_s = e*exp(s_s)*(n/n0)^(gamma_s-1)/(gamma_s-1),
    p_s = (gamma_s-1)*n*h_s,
    grad(p_s)/n = gamma_s*grad(h_s) - h_s*grad(s_s).

Here `s` is the dimensionless logarithmic entropy variable defined by the
existing EOS; `h` is specific internal energy, not specific enthalpy. Its
enthalpy is gamma*h. At fixed density the drift rate is
`C_s,t=-a_s*J dot grad(s_s)/e`. Grid thermal work uses the actual derivative
`dH/dC_s=V*h_s`, retaining both species' distinct local temperatures.

## A common physical projection

S samples regular fluid scalar coefficients at common physical quadrature;
N samples compatible field potentials there; W is cylindrical volume weight;
V is the existing fluid nodal volume. Define

    T = (N^T W N)^(-1) N^T W S,    D = P G T,
    A(J) = S^T W diag(J_h dot) D,    K(J) = (A-A^T)/2.

G is the compatible incidence gradient and P samples electric coefficients.
Tensor-product Gram solves implement T. This is a physical L2 projection of
the regular point entropy/enthalpy interpolants; it does not replace density or
alter the thermodynamic state. In general `S(C/n)` differs from `(S C)/(S n)`;
the latter is still used for independent continuous thermodynamic diagnostics.

The proposed point-sampling T in the earlier handoff gives a transpose with
large nonconverging local errors after skew symmetrization. Its unskew version
also has growing frozen modes. Both are retained as measured negative controls.
The physical projection supplies a consistent weak transpose in smooth cases.

## Skew current flux and pressure adjoint

The new pair is

    C_s,t = -a_s K(J) s_s / (e V),
    F_s = a_s/e { M_E G T(gamma_s*h_s - s_s*h_s/2)
              - P^T W [(S h_s) D s_s - (S s_s) D h_s]/2 }.

F_s is a weak electric forcing. For an Ampere current,
`J^T M_E G=0`. T and S preserve constants, so `A*1=0`, `1^T*A=0` and
`K*1=0` in exact arithmetic. Consequently each species conserves `sum V*C_s`,
constant entropy remains constant, and its grid thermal power equals `J^T F_s`.
Adding the inertial induction response gives zero grid thermal+current+magnetic
power. These are algebraic identities, not claims of local physical accuracy.

The graph charge flux `2*K_ij` is antisymmetric and has zero row divergence.
With species charge factor a_s/e, its central entropy flux is
`a_s*K_ij*(s_i+s_j)/e`. Thus each graph pair exchanges equal and opposite
entropy, including fluid elements with discontinuous traces. T is globally
supported through inverse Gram matrices and periodic Fourier interpolation;
these graph fluxes are not asserted to be exact current integrals through
geometric nearest-neighbor faces. A local physical enthalpy/Poynting interface
budget remains a coupled-adoption requirement.

For fixed current and density, K is skew and
`Q_s=sum V*n*s_s^2/2` is invariant. This gives a bounded finite-dimensional
linear evolution. It does not prove preservation of all entropy Casimirs,
nonlinear internal energy, positivity in continuous reconstruction, shock
entropy production or stability with evolving density/current.

## Analytic curl-kernel separation

Choose algebraic anchors ha=h[0], sa=s[0]. Write dh=h-ha and ds=s-sa. The
reduced electric term uses `-[(S dh)D ds-(S ds)D dh]/2`; the potential becomes
`T(gamma*h-s*h/2-ha*s/2+sa*h/2)`. This is the same equation. It isolates exact
compatible gradients before induction, avoiding cancellation of a large
pressure potential. The full electric field and original Ohm residual are
reconstructed afterward. Neither anchor is an equilibrium fit or thermal repair.

## Boundaries, time and approximation

The bounded cylinder has impermeable fluid walls, fixed axis/wall poloidal flux,
natural Btheta(R)=0, free wall Ez, and periodic z. There are no exterior plasma
cells, but this component supplies no evolving free plasma/vacuum interface.
Compatible weak charge continuity does not make the reconstructed current's
strong divergence, internal normal jumps or wall normal trace identically zero.
Independent local/axis/wall/max errors and strong divergence are retained.

Frozen-current transport uses explicit RK4 with finite positive EOS checks at
every stage; rejection leaves the input unchanged. It has no added numerical
diffusion. Dispersion, physical reconstruction error and finite-step RK4 drift
are separately measured. No shock limiter or entropy-production model is added.

## Background for the numerical construction

Skew forms, exact differentiation of constants, and compatible quadrature are
standard ingredients of entropy-stable discretizations; see Jesse Chan,
[Skew-symmetric entropy stable modal discontinuous Galerkin formulations](https://arxiv.org/abs/1902.01828).
The particular transfer, thermodynamic adjoint and limitations above are derived
for this repository and verified independently; that paper is not a proof of
this plasma scheme or of its coupled nonlinear stability.
