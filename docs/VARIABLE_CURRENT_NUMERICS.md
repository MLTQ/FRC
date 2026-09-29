# Variable-density compatible current and pressure-work bridge

This is a bounded instantaneous component toward the moving finite-mass solver.
It does not replace `SmoothFiniteMassMHD`, evolve a plasma-vacuum interface or
establish formation. Physical Ar+ and electron masses, quasineutral density,
axisymmetric periodic z, fixed radial magnetic flux, zero wall Etheta and natural
Btheta(R)=0 retain the preceding compatible-field contract. No diffusion,
resistivity, numerical atmosphere, floor or state-energy reset is introduced.

## One conservative density and the current-energy derivative

The caller owns n_cell=rho_cell/(mi+me), actual annular cell averages. Let R be the
degree-3/5 moment reconstruction and n=R n_cell at common positive physical
quadrature. The map preserves cell means, and the same **linear** R samples the
actual continuity rate n_t. It does not evolve another density. Negative means or
quadrature samples reject; positivity between samples and shock monotonicity are
not claimed.

P samples compatible electric/current cochains and W is cylindrical physical
quadrature. M=P.T W P and magnetic H are unchanged. Physical reduced mass
mu=mi me/(mi+me) gives

    I(n)=P.T W [mu/(n e²)] P,
    M J=−C.T H deltaB,
    K_J=J.T I J/2.

I is applied with tensor products and is symmetric positive definite for positive
n. The existing magnetic-space solver still checks the **original** mixed system

    M E=I J_t+F_other+F_n,
    M J_t=−C.T H B_t,
    B_t=C E,
    F_n=−P.T W [mu J n_t/(n² e²)].

Here J within quadrature expressions means P J. Current kinetic energy changes
with both current and density:

    K_J,t(cell)=integral_cell mu/(n e²)
               [J dot J_t−|J|² n_t/(2n)] dV / V_cell.

Its global density derivative in the actual cell-average coordinates is

    dK_J/dn_cell=R.T W[−mu |J|²/(2 n² e²)].

This adjoint matters: substituting point values or a nodal average would change
the work against the supplied mass flux. The new implementation exposes current
and density contributions separately. Full local B dot B_t includes the uniform
bias cross term, even though its volume integral vanishes at fixed-flux walls.

Multiplying the original mixed equations by current gives

    d(W_B+K_J)/dt=−J.T F_other
                 +integral mu |J|² n_t/(2 n² e²) dV.

The residual density term must exchange with bulk energy/pressure in the full
moving system. It is neither an extra physical heater nor permission to reset
conserved energy. This component does not yet supply that complete moving update.

## Compatible charge and the deliberate pressure bridge

The old dual graph's charge stream used a diagonal magnetic mass. With the new
physical H_g, its consistent replacement is S=(H_g g)/(dx dz), giving radial
flux I_r=S_previous_z−S and axial flux I_z=delta_r S, with zero endpoint streams.
These are weak graph fluxes; they are not claimed exact physical dual-face
integrals. They satisfy I_r=(M J)_r/dx, I_z=(M J)_z/dz and exact graph divergence
cancellation. Therefore their work against the compatible scalar gradient is
identical to physical magnetic/current work.

`CompatibleSpeciesPressure` explicitly keeps the older lumped nodal thermal
energy and reciprocal-logarithmic entropy mean:

    U=n w, a=(n/nref)^(gamma−1), A=a/(gamma−1), p=a U,
    W_s=sum V_dual A U,
    U_t=div_graph(w_face I/e),
    E_s=−grad(p/n)/e−w_face grad A/e.

It separates a complete scalar potential before the induction solve and returns
the full E for work. The closed fixed-n graph obeys W_s,t=J.T M E_s, conservation
of U, and zero logarithmic entropy rate. Physical electrons multiply this pair
by beta=mi/(mi+me), ions by −alpha=−me/(mi+me). Each species' work is verified
separately so a small mass fraction cannot hide an error.

**This pressure bridge is not high-order thermal closure.** Its nodal n is not
silently interchanged with R n_cell, its energy remains lumped, and nonlinear
face products retain pressure errors. Using the old stream with the new current
also gives an independently measured work mismatch. The bridge isolates those
two issues: correct common work is necessary but not sufficient for physical
accuracy. Full density compression, thermal energy placement, bulk advection,
vorticity, stress and complete current/enthalpy/Poynting flux remain to couple.

## Verification boundary

Focused tests include analytic polynomial density/current/field solutions,
conservative means and adjoints, central directional energy differences,
independent original three-block solves and positivity rejection. The smooth
nonpolynomial study supplies n_t from **surface-integrated prescribed particle
fluxes**, checked independently against continuum divergence. Its source is an
explicit manufactured instantaneous induction forcing, not a claimed parent-
species solution. Independent root-product field and primitive-polynomial
density sampling use increased quadrature to audit local energy and its rates.
Full parent-species and time-trajectory checks remain mandatory before adoption.
