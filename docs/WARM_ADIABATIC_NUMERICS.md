# Moving ideal electron and ion temperatures

## Scope and governing equations

This filled, bounded verification model assumes singly charged ions,
massless quasineutral electrons, scalar pressures and **common velocity**.
It implements the smooth ideal limit of the species thermal equations in the
[NRL Plasma Formulary](https://www.govinfo.gov/content/pkg/GOVPUB-D210-PURL-LPS118567/pdf/GOVPUB-D210-PURL-LPS118567.pdf),
printed pp.36–39, with zero heat flux, collision exchange and resistivity.
Both species have the same specified gamma, usually 5/3:

    rho = mi n,  v = vi = ve
    ∂t n + div(n v) = 0
    ∂t ps + v·grad ps + gamma ps div v = 0, s=i,e
    p = pi+pe, epsilon = p/(gamma-1)
    Etotal = rho |v|²/2 + epsilon + |B|²/(2 mu0).

Momentum, induction and total energy use the existing conservative ideal MHD
equations. With common velocity, summing species energy gives the existing
total flux

    Ftotal = [rho |v|²/2 + gamma(pi+pe)/(gamma-1)] v
             + (-v cross B) cross B / mu0.

There is no extra thermal-energy source to add: electron enthalpy is already
present once. This is an ideal common-velocity gate toward the full
two-temperature Hall system, not a replacement for its current enthalpy flux,
pressure battery, collisional partition or finite-mass terms. The exact
continuum tests have uniform axial B and fields depending only on z, so J=0
and curl(grad pe/n)=0. Species transport can therefore be tested without
conflating the missing Hall/pressure induction with the reference solution.
Their full pressure electric field generally has nonzero Ez along the radial
wall. The present common-velocity solver imposes only the **ideal MHD** electric
trace, so these fixtures are not a verification of a full conducting-wall
two-fluid electric boundary. That distinction must be resolved by the next
warm Hall boundary choice even though the omitted gradient is curl free here.

## Conservative mass placement and proof

Use x=r²/2, dV=2pi dx dz. Primal cell mass is Mij=Vij nij. Each of its four
corner dual volumes receives Mij/4; periodic z endpoints are folded. This
positive overlap map Q includes half supports at axis and wall. Electron
number density is always Q(rho)/mi, not a second independent field.

Let Fr,Fz be integrated particle fluxes from the **actual fluid Rusanov
step**, in particles/s. In each primal rectangle construct the RT0 flux:
its x component interpolates linearly between radial normal fluxes and is
constant in z; its z component is linear in z and constant in x. It matches
each normal flux and its divergence is constant within that rectangle.
The integral over a quarter rectangle is one quarter of the primal divergence.
Summing incident quarters cancels internal fluxes and proves

    Q(div_primal F) = div_dual(MapFlux F).

Interior radial dual flux is one quarter of the four surrounding radial
primal fluxes. At axis/wall average the two adjacent axial half faces.
Axial dual flux is one quarter of the surrounding axial primal fluxes; omit
absent exterior radial cells without renormalizing. This identity includes
general radial boundary flux. The warm fluid solver requires zero mass flux
at its reflecting slip radial wall and regular axis; z is periodic.

P, the Q1 corner-to-cell scalar average, is the volume adjoint of Q:
sum(Vp q P f)=sum(Vd Q q f). Overlap density is not an exact point sample.
Actual dual-volume centroids, rather than geometric corner locations, are
used for bulk reconstruction; their difference matters near boundaries.

## Entropy variable and compression

Define U=pe/(n/nref)^(gamma-1), w=U/n. Smooth ideal species equations imply

    ∂t U + div(U v)=0, Dw/Dt=0,
    pe=(n/nref)^(gamma-1) U,
    pe_dot=(n/nref)^(gamma-1) Udot+(gamma-1) pe ndot/n.

The second term gives compression work. For spatially uniform w, pe is
proportional to n^gamma, including discretely because U and n use the same
particle flux at each RK stage. A reference-density change rescales U only.

`bulk_entropy` uses Gamma times upwind reconstructed w. Order 1 is donor cell;
order 2 uses MC slopes in x,z and positive face-trace scaling. This scalar
scheme shares the fluid method's numerical mass diffusion. SSPRK2 advances
fluid and U with the same two stage states and mass fluxes. Both trial species
partitions must be positive before commit. There is no floor, adaptive retry,
temperature cap or total-energy reset. A fluid CFL bound is checked; it is
not claimed to ensure species positivity for arbitrary data.

## Species energy, entropy and numerical error

Cell electron internal energy is P(pe)/(gamma-1). Ion pressure is recovered
as conservative total thermal pressure minus P(pe). The adjoint identity
equates global dual and primal electron energy. Equal gamma keeps total fluid
dynamics independent of partition; the new path reproduces the existing ideal
fluid state bit for bit in tests.

Finite-grid mixing and total-energy dissipation do not satisfy an exact separate
ion entropy equation: residual thermal energy goes to the recovered ion
component. This is a numerical partition error that must converge in smooth
tests, **not** a model of shock or collisional heating. The linear integral
sum(Vd U) is conserved, but physical electron entropy is

    Se = kB/(gamma-1) sum(Vd n log[U/(n e)]) + constant * total particles.

The reconstruction/RK scheme does not guarantee a discrete physical entropy
inequality. The main strong contact has a small negative entropy change; its
refinement controls are retained. Exact total energy and positive temperature
alone do not establish entropy stability. Collision exchange, resistive
species heating, conduction and shock robustness remain missing.

## Interfaces and next integration

`ToroidalMHD._rhs_with_fluxes` exposes integrated mass flux without changing
the `_rhs`/`rhs` contract or fluid arithmetic. `DualMass` handles overlap and
flux integration, `bulk_entropy` only transport, `SpeciesThermodynamics` only
pressure/energy/diagnostics, and `AdiabaticSpeciesMHD` the coupled stages.

Next couple this moving state to the verified fixed-density electron-current
operator using a **local conservative** enthalpy/Poynting flux and mutually
consistent magnetic boundary trace. Do not reset total energy by the electron/
magnetic change. The thermal wall trace and conducting Hall wall are not
interchangeable. Full warm Hall waves and nonlinear evolution are still
required, followed by physical heating and finite-mass current transport.
Real vacuum, gas ends, heated edges, spatial coils and formation/UI goals
remain open.
