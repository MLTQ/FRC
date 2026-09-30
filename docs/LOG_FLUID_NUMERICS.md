# Conservative motion with physical log entropies

This bounded, smooth hydrodynamic component advances motion self-consistently.
It is a step toward replacing the older inaccurate moving pressure partition,
not a complete plasma solver or formation experiment. It retains both species
masses, M=mi+me, actual number density and separate physical pressures.

## State and parent equations

Stored cell means are n, m_r=n v_r, m_z=n v_z and L=n r vtheta. Two continuous
nodal scalars are s_a=log[p_a/(n e (n/n_ref)^(gamma_a−1))], n_ref=1e20 m^-3.
Thus p_a=n e exp(s_a)(n/n_ref)^(gamma_a−1), and H_a=p_a/(gamma_a−1).
The positive exponential is a physical variable transformation, not a floor.

The parent equations in an axisymmetric closed cylinder are

    n_t + div(n v) = 0
    (M n v_r)_t + (1/r) d_r[r(M n v_r²+p)] + d_z(M n v_r v_z)
        = (p+M n vtheta²)/r
    (M n v_z)_t + (1/r) d_r(r M n v_r v_z) + d_z(M n v_z²+p) = 0
    (M L)_t + (1/r) d_r(r M L v_r) + d_z(M L v_z) = 0
    (s_a)_t + v dot grad(s_a) = 0
    E_t + div[(K + sum_a gamma_a H_a) v] = 0,
    K = M n |v|²/2, E=K+H_e+H_i, p=p_e+p_i.

The last equation defines a separate conservative energy account. The evolved
physical state comes from momentum and the two entropies; temperatures are never
replaced by a partition of the account. The discrepancy is measured locally
and globally and is not claimed to vanish identically.

## Reconstruction, derivatives and boundaries

The existing density reconstruction preserves actual annular cell means.
Physical mr and L use regular r P(r²/2) and r² P(r²/2) polynomials; mz uses the
same map as n. Physical angular momentum divided by r gives n vtheta. Strong
interior velocity derivatives differentiate these same moment polynomials and
the quotient by n; cylindrical curl/divergence include the regular metric terms.
These point derivatives do not replace the conservative continuity law or
represent distributional jumps between polynomial cells.

Momentum and energy face fluxes use centered degree-3/5 normal moment traces,
mean-preserving transverse reconstruction and physical quadrature. Pressure is
evaluated from the face n and continuous entropy, not an independently fitted
pressure field. Axial/angular momentum and the energy account telescope.
The radial hoop force is integrated over each annulus. Axial boundaries are
periodic; axis and radial wall have zero normal particle flux. Radial pressure
traction remains. This is a closed fluid fixture, with no exterior matter.

One existing CompatibleParticleFlux drives continuity and both physical thermal
advection operators. Its integrated particle flux agrees with the mechanical
face quadrature. True physical K and both H are differentiated with the actual
n_t/m_t/s_t. RK4 also integrates their difference from conservative energy flux,
separating spatial defects from finite-step chain-rule error. All stages reject
nonfinite or nonpositive reconstruction atomically, without retry or repair.
There is no added numerical viscosity/diffusion, filtering, density atmosphere,
shock limiter, positivity theorem or entropy-stability theorem.

## Why a force-only energy correction is insufficient

For m_z=U n_cell, physical v_z=U is uniform under the same reconstruction.
Physical K gradients are dK/dm_z=M U times cell volume and dK/dn=−M U²/2
times cell volume. Any update conserving total m_z and n has integral K_t=0.
The previous bulk thermal transport nevertheless gives a small, nonzero
integral H_t for a general translating density/entropy field. Therefore an
exact energy-canceling pressure adjoint would require a spurious net force.
Fixing pressure alone cannot guarantee both energy and momentum conservation.

The independent translation probe measures this defect versus mesh, flow sign
and quadrature. The conservative stress is retained and physical-energy error
is exposed in trajectories. A joint transport change remains necessary if its
convergent error is unacceptable; the energy account is never used as a repair.

## Evidence and remaining scope

Independent exact solutions include a translating two-temperature contact,
a finite-amplitude isentropic simple wave before shock formation, and solid
rotation with dp/dr=M n Omega²r. Tests check their original Euler/pressure
equations, actual unforced evolution, regular-axis curl and physical kinetic
work. The study retains initial projection, all accepted histories and five
snapshots, with separate higher-order physical audits and bitwise replay.

Field/current evolution, Maxwell/drift stress, current enthalpy/Poynting work,
collisions, physical transport, heated free edge/exact vacuum, real ends/coils,
shock/kinetic validity and converged FRC formation remain unfinished. Neither
these fluid fixtures nor the previous initialized FRC viewer establishes
three-dimensional stability or particle confinement.
