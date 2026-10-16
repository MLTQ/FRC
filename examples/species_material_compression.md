# Smooth neutral compression and rebound

`run(destination)` evolves neutral argon at n0=1e19/m³ in |X|<a=25mm,
T0(X)=1eV[1-(X/a)²], with an initially inward affine velocity. The physical
Euler solution has x=s(t)X, n=n0/s, v=s_dot X and T=T0(X)s^(1-gamma).
The pressure gradient gives s_ddot=A s^-gamma, A=2eT0/(ma²), gamma=5/3.

The initial velocity is chosen so the continuum solution turns at s_min=.5:
s_dot(0)=-sqrt(2A/(gamma-1)[s_min^(1-gamma)-1]). The independent first
integral is s_dot²=2A/(gamma-1)[s_min^(1-gamma)-s^(1-gamma)]. Scalar quadrature
with s=s_min+y² removes the square-root turning singularity and supplies
the full compression/rebound duration. At that duration s=1, the velocity
is the opposite of its initial value, and temperature returns to its initial
profile. The numerical dynamics are never constrained to this trajectory.

Runs all32/64/128 cosine meshes with C1=.5,C2=2, comparing the original
compression stress and the minmod-limited stress. This specifically checks
whether unnecessary numerical heating in smooth compression is removed while
the complete solution and energy remain accurate. It is an independent
continuum reference for the concrete viscosity change, not a unit test or
generic verification framework. No charge, field, collisions or reactions.

Every accepted state, particle amount, energy history, source identity and
failure are saved. Endpoint errors refer to the independent full-cycle
reference and are meaningful only for successfully completed cases. Exact
vacuum, positive actual thermal energies and original physical mass are retained.
No failed history is shortened, retried with hidden settings or repaired.

Usage: `python -m examples.species_material_compression FRESH_OUTPUT` with
the usual headless plotting environment. PNG/PDF show the physical front
against independent quadrature and the actual conservation history.

`--shock` instead runs a finite free neutral argon slab with the existing
`mhd_benchmarks.sod_exact(..., gamma=5/3)` central Riemann solution. Initial
left/right number densities are1e19/1.25e18m^-3 and temperatures1/.8eV;
a=50mm, speed0=sqrt(eV/m), duration=.12a/speed0. Both exterior regions are
exactly empty, with evolved material fronts. Only |x|<a/2 is compared, before
edge rarefactions reach that region. The whole-domain mass, momentum and energy
remain physical; there are no clamped support velocities or ghost particles.

All64/128/256 meshes run with both original and limited stress. Central
normalized density/velocity/pressure L1 differences use16-point Gauss
quadrature split at physical cells and analytic fan/contact/shock boundaries.
This shock comparison complements the smooth rebound: lowering viscosity must
retain a credible shock solution. Raw histories, numerical heating and actual
energy/momentum are saved; PNG/PDF show the finest actual fields against the
existing independent Euler solution. No unit tests are invoked.
