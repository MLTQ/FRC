# Warm expanding gas with actual neutral radiative transfer

This is an independent coupled motion/source benchmark for common density
coarsening. A declared mathematical two-level atom has hydrogen mass1.00794u,
excitation3eV and spontaneous decay rate2e5/s. Both levels are neutral and have
the same mass. Initial excited gas has n=1e19/m³,half-width a=25mm,velocity0,
and T=.1eV[1-(X/a)²]. Fine-cell thermal values are exact quadratic averages.
There are no electrons, ionization/recombination or electric field in this
specific soluble physical limit; exact exterior vacuum remains empty.

Decay transfers N, momentum and energy to the lower level and emits3eV photons
per transition, leaving its donor velocity and translational temperature intact.
Both fluids therefore share the independent affine Euler solution x=sX,
v=s_dot X,n_total=n0/s,T=T0 s^-2/3[1-(X/a)²], with
s_ddot=A s^-5/3,A=2eT0/(m a²). The scalar energy-integral quadrature used in
atomic_moving.md with plasma2=0 fixes the FULL time to s=1.5 and its speed.
It does not integrate the production material equations to generate a reference.
N_upper=N0 exp(-2e5 t), photons=3eV(N0-N_upper), and total matter/binding/photon
energy is constant. Photon transport/recoil remain the source model assumptions.

All cases use source-half/moving-full/source-half with actual new lower-level
supports and populations, four initial density cuts per mechanical element,
physical mass, and disclosed minmod Q=(.5,2). The16/32/64 mechanical grids use
density width fraction.1; additional32-cell fractions.05/.025 examine geometric
coarsening sensitivity. Protected support endpoints are never merged.
No thermal repair, clipping, assigned product velocity or shortened success.

All accepted ragged states, normalized N/U scales, physical budgets, input/source
identity and PNG/PDF are retained. Resource bounds10000 cells/steps explicitly
fail the original duration. Failed states retain their last overlap and do not
receive final analytic endpoint comparisons. This benchmark does not replace
the still-required full hot charged argon attempts or any FRC formation gate.

Use `python -m examples.reactive_material_affine FRESH_OUTPUT`. No unit tests.

First retained case in outputs/reactive_material_affine_coarsened fails after
151 steps at4.158318µs of6.083704µs. The process subsequently disappeared without
a recoverable exit output; termination.json distinguishes that from its stale
summary flag. No missing cases are counted as successful. Replaying the last
saved state reproduces collapse at RK4 k4: distinct support-edge cells should
have width3.479932e-18m under the fixed affine map, but both endpoints round to
the same absolute coordinate after crossing0.03125m, where the ULP doubles to
6.938894e-18m. Coarse widths remain positive. geometry_diagnosis.json retains
hex coordinates, the input hash and independent Decimal width arithmetic.
This is a geometry-precision defect (frc-transport.5.5.3), not evidence for
adopting coarsening or physically collapsing the support cells.
