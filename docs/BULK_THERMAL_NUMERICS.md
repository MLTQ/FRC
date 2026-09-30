# One continuity flux for moving density and physical entropy

The physical log-entropy component now has bulk transport through the actual
conservative continuity flux. It advances changing density and both species
entropies in prescribed-flow verification. Momentum, magnetic/current dynamics
and complete plasma total energy are not yet advanced by this new component.
The older full moving solver retains its source and measured accuracy limits.

## Conservative cells and continuous normal particle flux

Use x=r²/2, with physical volume 2 pi dx dz. Cell means store n and the actual
particle momenta n vr and n vz. One shared reconstruction gives face means

    qx = r n vr,    qz = n vz,
    n_cell,t = -Dx qx - qz Dz.T.

Here Dx and Dz are conservative differences divided by the respective cell
widths. Axial faces are periodic; qx=0 at the regular axis and impermeable wall.
Multiplication by 2 pi dz or 2 pi dx gives actual particles/s through the faces.
The radial momentum fit uses its regular r P(x) form. No analytic velocity or
numerical current is fitted into these face values.

The normal face fit uses ONE stencil centered on the face and an integer face
index. The initial attempt reused the arbitrary-point density sampler's average
of two cell-containing stencils. At a normal face this biases advection and
causes positive spatial growth; full-transit density became nonpositive.
The centered correction retains the existing density reconstruction and flux
primitive. In periodic constant flow its operator is skew-symmetric, with no
spatial growth in either direction. No added diffusion masks the instability.

The existing density reconstruction R uses degree-p cell polynomials that
preserve their actual means. For one coordinate define the flux primitive

    (H q)(x) = q_left + integral_left^x R(D q)(a) da.

Mean preservation gives both actual face traces and continuity of the normal
flux. Differentiation gives dH=R D. Tensor products then supply

    Gamma_r = (Hx qx Rz.T)/r,
    Gamma_z = Rx qz Hz.T,
    div Gamma = -R(n_cell,t).

Density may have jumps between the reconstruction's local polynomials, but
normal particle flux is continuous. This is a reconstruction of one continuity
law, not a separate nodal density evolution. The independent verification
constructs H by root-product interpolation of cumulative face values instead
of the production moment solve and primitive integration.

## Bulk log entropy and its pressure work

With n=R n_cell, nodal s, scalar sampling N and the same physical metric
M_n=N.T W n N as the current-driven thermal component, define

    M_n s_t = -N.T W [Gamma dot grad(Ns)].

Constant entropy is preserved exactly. The physically relevant total entropy
is proportional to S=integral n Ns dV. Integration by parts, continuous normal
Gamma, continuous Ns and the shared continuity equation give

    S_t = integral [n N s_t + R n_cell,t Ns] dV = 0.

Both terms are retained when density changes. This is a spatial work identity
up to solve, quadrature and floating-point errors; it is not a fully discrete
entropy theorem. Gamma/n is the transport velocity of this lift. It need not
equal the pointwise ratio of independently reconstructed momentum and density
at finite resolution; that distinction is essential for future mechanical work.

Pressure and internal energy retain the existing physical convention

    p = n e exp(Ns) (n/nref)^(gamma-1),    epsilon = p/(gamma-1),
    epsilon_t = epsilon [N s_t + gamma R n_cell,t/n].

Let h solve M_n h=N.T W epsilon and let g_n=R.T W gamma epsilon/n. Write
n_cell,t=T q and Gamma=L q, where T includes the negative divergence. The
thermal pressure-work coefficient on admissible faces is

    w_p = T.T g_n - L.T W [(N h) grad(Ns)],
    dH_thermal/dt = q dot w_p.

`face_work` is this generalized coefficient in the explicitly stated qx/qz
units. It is not a momentum force or an unexplained heating term. In the
continuum its global power equals -integral p div v dV at the closed boundary.
No thermal energy is fitted to a chosen force or reset after a step.

The existing species current-driven rates may be added with beta/-alpha factors;
they retain their separate compatible pressure-field work. This study isolates
bulk transport and compression and does not yet claim the full nonlinear sum.

## Physical quadrature and numerical error

For degree p in x, the entropy integration-by-parts products reach degree 2p in
x, hence degree 4p+1 in r after the cylindrical volume factor. Physical radial
Gauss order must therefore be at least 2p+1. The new bulk component rejects
underresolved quadrature instead of claiming its weak identity. Default p=5
uses order 12, with higher independent controls. Exponential thermodynamics is
not polynomial, so energy quadrature convergence is still measured separately.

There is no added viscosity, Rusanov dissipation, limiter, filter, carrier floor
or numerical atmosphere. The unlimited smooth interpolation can overshoot and
has no shock-admissibility theorem. Reconstructed nonpositive density or
unrepresentable exponential pressure rejects without repair. Truncation,
projection, interpolation and finite solve precision cause the measured errors.

## Evolving-density verification and atomic time integration

The prescribed-flow driver advances n_cell, s_e and s_i together with RK4.
Each stage's particle momentum is computed from the current numerical density
and a declared velocity, never from the exact evolving reference density.
All stages and the endpoint must be finite with positive reconstructed density
and species pressure. A rejected step leaves input arrays intact.

The two reference flows are uniform translation for a complete axial transit,
and v_z=V sin(kz) for |k V t|=ln 2. The latter has inverse characteristics
tan(k xi/2)=exp(-k V t) tan(kz/2), n=n_initial(xi) dxi/dz and s=s_initial(xi).
It tests real compression and adiabatic temperature change. It is not a solution
of the unforced momentum equations; the velocity is explicitly prescribed.

A separate RK4 ledger integrates q dot w_p. Physical internal-energy change
and that ledger remain independently computed and their discrepancy is measured.
Since S is bilinear in evolving n and s, RK4 does not generally preserve it
exactly. Time controls must quantify its finite-step drift. Particle number is
a linear invariant of the telescoping continuity law.

## Remaining moving-state and mechanical contract

This work resolves the shared density/bulk-entropy transport portion of the
state contract. It does not justify dropping two nodal log entropies into the
old cell total-energy/ion-pressure subtraction. The next candidate state must
explicitly retain true conservative n and all momentum/angular-momentum means,
affine magnetic flux/current, and both physical thermal fields. A separately
evolved conservative cell total-energy account may diagnose local discrepancies;
its conservation alone cannot establish conservation of the physical components.
It must never be used to reset pressure or silently replace a component energy.

For reconstructed physical bulk kinetic energy K(n,m), its momentum derivative
is the volume pullback of velocity, not generally the cell-center m/rho. The
pressure part of the mechanical update must satisfy

    K_m dot m_t,pressure = -q dot w_p,

while retaining the correct conservative axial momentum, angular momentum and
cylindrical stress balance. Merely applying the face adjoint as an acceleration
does not prove those properties. The next implementation must derive this local
stress/work pairing with the true kinetic metric, velocity/vorticity and shared
mass flux, then verify independent polynomial forces, uniform flow, translating
two-temperature contact and all finite-mass wave branches. Full current/drift
stress and enthalpy/Poynting energy transport must then use that same contract.

The variational work-pairing approach is related to the broader fluid methods of
[Gawlik and Gay-Balmaz](https://arxiv.org/abs/2211.08745) and classical MHD methods
of [Gawlik, Gay-Balmaz and Manach-Pérennou](https://arxiv.org/abs/2503.01239).
Our specific mixed conservative-cell/log-entropy construction is derived and
tested here; their full-method conservation and thermodynamic results do not
automatically apply to it.

Heated free-edge/exact-vacuum coupling, physical collisions/transport/ionization,
real ends/spatial coils, shocks and kinetic validity, converged FRC formation and
the interactive bench remain unfinished. Axisymmetry does not prove 3-D stability
or particle confinement.
