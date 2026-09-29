# Conservative warm massless-electron Hall reference

`WarmHallMHD` couples moving ions, all axisymmetric magnetic components and
separate equal-gamma electron/ion pressures in a filled periodic cylinder.
This is an experimental bounded reduction. The pressure-only splitting has a
new high-resolution rejection and a continuum principal-symbol limitation;
it is not a generally validated warm formation solver. See
WARM_HALL_VERIFICATION and WARM_PRESSURE_SPLIT_LIMITATION.

No magnetic vacuum, physical gas ends, coils, collisions, conduction, ionization
or electron inertia are included. Resistivity is explicitly rejected pending
a physical heating partition. The original low-density argon formation target
is unchanged; this bounded verification fixture does not establish its validity.
Its base n=1e20 m^-3, R=.035 m, Te=5 eV and Ti=1 eV match the original input;
uniform equilibrium and idealized boundaries are the verification simplifications.

## Parent conservation

For quasineutral singly charged species, rho=mi n, v=vi, ve=v-J/(n e), and
J=curl B/mu0. With isotropic pressures, equal gamma and no physical dissipation:

    rho_dot + div(rho v) = 0
    (rho v)_dot + div[rho vv + (pi+pe+B²/(2mu0)) I - BB/mu0] = 0
    E = -v cross B + J cross B/(n e) - grad pe/(n e)
    B_dot = -curl E
    epsilon_e_dot + div(epsilon_e ve) + pe div ve = 0
    epsilon_i_dot + div(epsilon_i v) + pi div v = 0.

epsilon_s=ps/(gamma-1). For K=rho v²/2 and total energy
Et=K+epsilon_i+epsilon_e+B²/(2mu0), the flux is

    Ft = (K+epsilon_i+pi) v + (epsilon_e+pe) ve + E cross B/mu0.

Relative to the ideal common-velocity flux, the increments are

    FH = [J cross B/(n e)] cross B/mu0
    Fp = [-grad pe/(n e)] cross B/mu0
         - gamma/(gamma-1) (pe/n) J/e.

Both terms in Fp are retained, including a pressure electric gradient with zero
induction curl. Cell total energy is never reset to match a separately computed
electron/magnetic change. Momentum includes total pressure and Maxwell stress
once; no extra Hall force is added to ion momentum. The parent species equations
are in the [NRL formulary](https://www.govinfo.gov/content/pkg/GOVPUB-D210-PURL-LPS118567/pdf/GOVPUB-D210-PURL-LPS118567.pdf),
printed pp.36–39. Finite-mass current energy/transport is a subsequent extension.

## Unified bounded magnetic trace

At the radial wall psi is fixed, Br(R)=0, Etheta(R)=0 and **Btheta(R)=0 with
free Ez(R)**. Axial boundaries are periodic and the fluid wall is impermeable.
This is not the old conducting Ez=0 condition or a plasma–vacuum interface.
Wall pressure/density may vary and produce tangential Ez. Zero Btheta and
Etheta give zero radial Poynting power pointwise; constant zero Btheta also
gives zero normal charge/electron enthalpy flux.

`TraceHallCurl` extends Ez through the wall and rebuilds its positive mass
matrix and Faraday boundary column. In x=r²/2,z, Ez is Q1 in x and constant
in z, while Er/r is constant in x and Q1 in z. The full nodal Q1 scalar
gradient belongs to this electric space, and C grad=0 at axis, wall and seam.

Magnetic energy uses the compatible psi form plus sum V <r²> g²/(2mu0), with
g=Btheta/r constant per cell. The weak current is M j=-C.T H b. Its poloidal
work agrees with every DOF of the conservative dual current constructed from
K=2pi <r²>g/mu0 and actual Kwall=0. The independent zero-J1 toroidal diffusion
mode tests this boundary. Evolved fields converge, but the instantaneous
first-annulus resistive derivative remains inaccurate; it is not hidden by
global energy conservation or the passing evolved-mode test.

## Conservative electron-current enthalpy flux

`DualMass` overlaps actual ion mass into dual nodal volumes; no independent
electron density is evolved. The Hall denominator is its positive Q1 sample.
Set U=pe/(n/nref)^(gamma-1), with nref=1e20 m^-3 in this protocol. The ideal
bulk step transports U with the commuting dual mass flux. Pressure/current
steps use the means in ELECTRON_PRESSURE_NUMERICS. Electron cell energy is
the exact Q1 average of pe/(gamma-1), adjoint to the mass overlap.

Within each dual rectangle, reconstruct current with a radial component linear
in x/constant in z, and an axial component linear in z/constant in x. Its
normal face values match the integrated charge flux and its divergence is zero.
Integrating through primal faces preserves local zero divergence for constant
scalar weight, including axis/wall half supports.

Multiply by primal Q1 temperature Theta=pe/n [J]. Each primal face crosses
two half dual cells. Exact integration of the half linear temperature yields
3/4 and 1/4 endpoint weights. Their sum gives the actual current-temperature
face integral; multiply by -gamma/[(gamma-1)e] to obtain watts.
`dual_enthalpy` implements this integration, not an interpolated energy change.

`electron_energy` separately integrates pressure Poynting power. With
Ep=(r er,0,ez), Bt=r g, the physical poloidal flux is (-ez Bt, er r Bt)/mu0.
Radial traces use the actual wall Bt; axial averages retain the exact annular
<r²> moment. Ep crossed with poloidal B is toroidal and has zero axisymmetric
divergence. Summed face powers update Et by conservative volume divergence.
The Hall step uses its full compatible Poynting flux at the same boundary.

Global total energy is conserved, but the **local** flux divergence does not
exactly equal mapped electron-plus-magnetic energy change. Their difference
is numerical ion-energy redistribution, not collisional heat. Its convergence,
maximum local value and high-resolution rejection are explicitly measured.

## Time integration and scope

Stages are ideal(dt/2), pressure(dt/2), Hall(dt), pressure(dt/2), ideal(dt/2).
Pressure steps hold n/momentum fixed and use magnetic midpoint and the
logarithmic thermal mean. Hall induction uses compatible implicit midpoint
and the actual stage density. Ideal steps use shared fluid/species RK2.
The composition is palindromic and second-order consistent for smooth resolved
solutions; explicit RK2 does not make the full update exactly reversible.

All five stages must pass species positivity and nonlinear residual checks
before accepted state, time or accounting changes. Only the computational
Hall preconditioner may change on rejection. No floors, energy repairs or
hidden retries are used. Fluid CFL alone cannot guarantee whistler phase
accuracy, pressure convergence or stability; explicit caps are tested.

Independent uniform-background warm modes use cs²=gamma(pi0+pe0)/rho0.
Since div J1=0, div ve1=div vi1, and each pressure perturbation is its
background fraction of total acoustic pressure. The pressure electric gradient
has zero linear induction curl but remains in energy transport. Nonlinear
states and finite-amplitude controls are separate checks.

The ideal bulk entropy defect remains. No shock entropy inequality, general
positivity theorem, physical shock heating partition or uniform accuracy
across disparate wave speeds is claimed. The isolated pressure split admits
a nonhyperbolic frozen symbol for this nonbarotropic reference. Joint Hall/
pressure integration is the next derivation and validation gate; smaller steps
or a numerical heat source alone cannot justify this split. Stage 2 stays gated.
