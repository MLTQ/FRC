# Conservative moving finite-mass reference

`FiniteMassMHD` advances density, all bulk momentum components, total energy,
magnetic flux/toroidal field and electron entropy through one unsplit derivative.
It is a bounded smooth-reference solver with physical electron mass. The
pressure errors in FINITE_MASS_VERIFICATION prevent treating it as an accurate
formation model. No production Stage 2 gate is enabled.

## Equations and assumptions

The exact reduction and independent parent-species tests are in
MOVING_CURRENT_NUMERICS. With M=mi+me, rho=M n, alpha=me/M, beta=mi/M,
delta=beta-alpha, mu=mi me/M and w=J/(ne),

    vi=v+alpha w, ve=v-beta w, J=curl B/mu0,
    rho_t+div(rho v)=0,
    (rho v)_t+div[rho vv+mu n ww+(pi+pe+B²/(2mu0)) I-BB/mu0]=Fi+Fe,
    B_t=-curl E,
    ps_t+vs dot grad ps+gamma_s ps div vs=(gamma_s-1) Qs.

Full E retains Hall, both pressure gradients, changing-density electron inertia,
bulk advection, velocity gradients and mixed kinetic scalar potentials. It is
the same physical finite-mass Ohm law documented and verified in that reference;
the actual continuity derivative now supplies n_t, rather than analytic n_t.

Let Kb=rho v²/2, Kj=mu n w²/2 and Hs=gamma_s ps/(gamma_s-1). Then

    W=B²/(2mu0)+Kb+Kj+pi/(gamma_i-1)+pe/(gamma_e-1),
    F=(Kb+Kj+Hi+He)v+E cross B/mu0
        +w[mu n v dot w-delta Kj+alpha Hi-beta He],
    W_t+div F=vi dot Fi+ve dot Fe+Qi+Qe.

Quasineutrality, isotropic species pressures and no displacement current are
assumed. No physical viscosity, collisions, resistivity, heat conduction,
ionization, radiation, neutral coupling or kinetic stress is present. Explicit
manufactured force/heat is derived from the parent equations; unforced tests
have no callback. A consistent me=0 comparison is available, without fitting an
effective mass. Different gamma_e and gamma_i remain explicit.

## State and thermal partition

Primal annular finite volumes store `(rho, rho vr, rho vz, W, r rho vtheta)`.
Corner psi and cell g=Btheta/r determine the compatible magnetic field; dual
nodes store Ue=pe/(n/nref)^(gamma_e-1). The same quarter-volume overlap maps
rho/M to dual n. Electron energy uses the volume-adjoint Q1 pressure average.

Bulk energy uses vr=mr/rho, vz=mz/rho, Omega=L/(rho <r²>) and
Kb=rho(vr²+vz²+<r²>Omega²)/2. Radial subcell kinetic variation is not reconstructed
from the cell mean; this is a declared spatial approximation. Magnetic and
physical current energy use compatible physical quadrature, including local
1/n. Ion pressure is the residual `(gamma_i-1)(W-Kb-WB-Kj-epsilon_e)`.
Its dual overlap supplies the ion pressure in Ohm/current enthalpy terms.

The initial state sets this numerical energy convention once and separately
measures its difference from independent parent-species energy. Evolution never
resets total energy to fit a pressure, nor adds missing thermal energy. Component
accounting is an algebraic consistency check, not independent proof of accurate
ion temperature. Rate and trajectory pressure oracles supply that separate test.

## Shared spatial derivative

1. Axis-regular limited second-order bulk reconstruction yields shared Rusanov
   fluxes of all five conservative variables. Separate electron fraction and
   current-energy traces account for unequal gammas and Kj. Maxwell stress is
   included once. The central energy flux contains bulk enthalpy and kinetic
   advection; computed Poynting/current powers are added separately.
2. The actual integrated mass flux gives primal rho_t, dual n_t and dual particle
   transport. Its overlap commutes with divergence. No second continuity law is
   numerically advanced. A mass-weighted regular velocity map supplies vr/r,
   Omega and vz for induction; the impermeable wall trace is explicit.
3. `MovingCurrent` solves the complete instantaneous finite-inertia problem.
   Its original mixed residuals retain tolerance 1e-10; explicit linear iterative
   refinement is permitted without changing that tolerance or the physical step.
4. Shared face powers add full E cross B/mu0 and mixed/current enthalpy transport.
   Drift stress transports mu w/e with the divergence-free dual charge flux;
   angular stress includes r wtheta and radial momentum has its hoop source.
5. Electron U uses the same bulk particle flux plus beta-scaled electron-current
   entropy transport and any explicit heat source. Inferred pi is checked at
   every stage. Local adiabatic ion consistency is measured, not exact by fiat.

Rusanov diffusion uses the fast MHD speed with gamma p=gamma_e pe+gamma_i pi.
It is **numerical dissipation**, not a physical heating partition or a complete
bound on dispersive two-fluid eigenfrequencies. Induction has no added diffusion.
Slope admissibility changes reconstructed traces, never cell averages. There
is no proven entropy inequality or full-update positivity guarantee.

## Time integration, boundaries and failure

Classical explicit RK4 evaluates this full derivative at all four stages.
An explicit caller-specified maximum time step is required. It is not a general
stability theorem; spatial/time controls and physical mode frequencies must be
resolved. All stages and the final state must have finite fields and positive
density/pe/pi before accepted state, time or integrated source/boundary work
changes. Rejected steps raise atomically, without hidden retries or floors.

The cylinder has periodic z, fixed psi on the axis/wall, Btheta(R)=0/free Ez(R)
and zero normal particle/current flow. The axis is regular. These reference
conditions exclude real gas ends and the plasma-vacuum interface. No atmosphere
is introduced: this verification domain is intentionally filled with plasma,
and eventual formation must preserve exact exterior vacuum.

## Diagnosed accuracy limit and next correction

Slow nearly solenoidal modes require radial and axial compression to cancel to
high precision. Limited central normal-momentum reconstruction spoils this
cancellation; the measured diffusive contribution is much smaller. Decreasing
the time step cannot remove a spatial mass/enthalpy error. Exact global total
energy also cannot prevent spurious local transfer into inferred ion pressure.

`annular_polynomial` separately fits actual cell moments of r^s P3(r²/2), with
regular powers s=0/1/2 and periodic axial cubic interpolation. Its continuity
diagnostic improves the slow-mode rate substantially, but it is **not used** by
FiniteMassMHD. It has no limiter or positivity theorem. The next work must
reconstruct shared particle, momentum, enthalpy and kinetic transport consistently;
changing only mass flux is not a valid thermal fix. Repeat local pressure-rate,
uniform-flow/contact, wave and independent conservation tests before adoption.

Physical closures, shocks and entropy production, kinetic ordering, heated free
edge/exact vacuum, real ends/spatial coils and converged formation remain open.
Closed flux or an axisymmetric equilibrium does not establish three-dimensional
stability or particle confinement.
