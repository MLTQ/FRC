# Reacting species and local electromagnetic feedback

`AtomicElectromagnetic` couples every actual heavy-level/electron momentum to
the laboratory electric field, with physical charge-state masses. This is the
local source of the nonrelativistic multifluid–Maxwell equations:

    m_s d v_s/dt = q_s (E + v_s cross B)
    epsilon0 d E/dt = -sum_s q_s n_s v_s.

Within this source step, populations and B are fixed; species thermal energy is
constant. Sources and spatial fluxes change them in their own equations. The
method follows the time-centred source approach described by
[Wang et al. (2020)](https://arxiv.org/abs/1909.04125); its implementation here
uses a three-by-three Schur solve rather than a full species matrix.

Let h_s=sqrt(m_s f_s/m0) v_s/v0 and e=sqrt(epsilon0/(Nref eV)) E, where
v0=sqrt(eV/m0). Then dh_s/dt=omega_s e+h_s cross Omega_s,
de/dt=-sum omega_s h_s, with signed
omega_s=q_s sqrt(Nref f_s/(epsilon0 m_s)), Omega_s=q_s B/m_s.
This block is skew-symmetric. Implicit midpoint preserves its quadratic energy
up to arithmetic/linear-solve error. With a=dt*Omega/2 and C h=h cross a,
R=(I-C)^-1=(I+C+a a.T)/(1+|a|²). Eliminate the midpoint h values to solve

    [I+(dt/2)^2 sum omega_s² R_s] e_mid
      = e_old-(dt/2) sum omega_s R_s h_old.

Recover each h_mid, then h_new=2h_mid-h_old and e_new=2e_mid-e_old. Populations
remain exactly unchanged. All empty species stay empty; neutral matter feels
no Lorentz force. No floor, artificial mass, renormalization, charge projection
or energy correction is used. Products are ordered with sqrt(f) so tiny physical
populations do not first underflow when multiplied by a kilogram mass.

The atomic state stores relative momentum P' and energy E' in a fixed frame U.
The field solve uses lab velocities. Its exact midpoint work transforms back as
Delta E'_s=h_mid dot Delta h_s-(U/v0) dot Delta P'_s. This integrates the actual
species work equation; it does not reset a thermal state to a chosen budget.
Laboratory current uses all actual charge numbers, including charge residuals.

`step` composes material-source half, electromagnetic full, material-source half.
The material source can include the existing reactions and charged collisions.
Heat samples the real quarter/three-quarter times; no populations or moments
are reset between operators. It is second order for smooth resolved sources.

`observe` includes electric energy and field momentum epsilon0 E cross B.
For a homogeneous charge-neutral system with constant B and no gradients,
particle plus electric energy (and binding/photons/input heat), and particle
plus field momentum, are independent conservation budgets. Constant magnetic
energy can be omitted from changes. An external current or varying prescribed
B would require additional source work and is not implemented here.

This includes displacement current and plasma oscillations, unlike the moving
solver's quasineutral magnetic/current constraint. Large dt is stable for this
linear source but has phase error; omega_p*dt and cyclotron*dt are reported.
It must NOT be spliced into the old constrained state or applied to arbitrary
spatial charge transport without compatible Maxwell curls and Gauss continuity.
No spatial field evolution, free-boundary change, wall, sheath or formation
claim follows from this homogeneous coupling. Moving assembly remains required.

Optional `method="exponential"` uses the same physical source matrix via the
sparse exponential action in species_field_exponential.md. This avoids midpoint
phase dispersion inside that substep. Species work is the endpoint kinetic
difference, [(h_new+h_old)/2] dot (h_new-h_old), transformed to the stored frame;
E is independently evolved, never adjusted to the particle budget. No global
energy/temperature repair occurs. Reaction/collision splitting accuracy still
requires full histories and timestep comparisons. Default midpoint is retained
as the independently measured earlier method.
