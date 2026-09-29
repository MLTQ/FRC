# Next bounded gate: separate electron and ion thermal transport

The compatible Hall solver now advances nonlinear compressible fields and flow,
but it uses pe=0 and a single represented gas energy. Its prescribed-resistance
heat is a conditional reduced closure, not the original Te=5 eV / Ti=1 eV
argon target. The next step must derive and verify massless electron/ion thermal
transport before extending the moving free surface or claiming formation.

Current progress: the conservative dual current and isolated electron pressure/
entropy operator are implemented and verified in ELECTRON_PRESSURE_NUMERICS and
ELECTRON_PRESSURE_VERIFICATION: 266 tests, 46 records, actual heat transport and
nonlinear pressure battery. The moving ideal common-velocity limit is now
implemented in `AdiabaticSpeciesMHD`: overlap mass/flux, shared fluid/electron
RK stages, compression and positive species energy partition. See
WARM_ADIABATIC_NUMERICS and WARM_ADIABATIC_VERIFICATION. The next integration is
**local electron-current enthalpy/Poynting transport with compatible warm Hall
boundary traces**, then physical heating/exchange and finite-mass current terms.

## Parent conservation to preserve

In the quasineutral, singly charged, massless-electron reduction, v=vi,
ve=v-J/(ne e), rho=mi ne. With scalar isotropic species pressures and no
viscosity for the first bounded reference:

    E = -v cross B + eta J + J cross B/(ne e) - grad pe/(ne e)
    ∂t epsilon_i + div(epsilon_i v) + pi div v = Qi - div qi
    ∂t epsilon_e + div(epsilon_e ve) + pe div ve = Qe - div qe.

Collision exchange must cancel between Qi and Qe. Their sum supplies positive
eta J² under the chosen collisional Ohm reduction; its species partition needs
a physical model, not an unexplained split constant. Additional modeled heating
or losses must have explicit external reservoirs and units. First verify an
ideal, adiabatic two-temperature limit independently; only then introduce
collisional exchange, resistivity and justified conduction/ionization.

With K=rho v²/2 and U=K+epsilon_i+epsilon_e+B²/(2mu0), the total energy flux is

    F = (K+epsilon_i+pi) v + (epsilon_e+pe) ve + E cross B/mu0
        + qi + qe.

Relative to a single common-velocity gas flux it includes
**-(epsilon_e+pe) J/(ne e)**. Adding the pressure term to Ohm and its Poynting
flux while omitting this electron enthalpy drift is not a valid thermal closure.
For variable density, curl(grad pe/ne) generally does not vanish. A barotropic
pressure whose induction curl vanishes still carries electron pressure work
and enthalpy. Momentum uses total pi+pe and the existing Maxwell stress once.

These equations are a derivation target under the stated reduction, not an
already implemented or validated physical regime. Parent species equations
and transport assumptions are in the [NRL formulary](https://www.govinfo.gov/content/pkg/GOVPUB-D210-PURL-LPS118567/pdf/GOVPUB-D210-PURL-LPS118567.pdf),
printed pp.36–39; finite-mass terms are separately required by
[Martin's generalized current-energy derivation](https://ecommons.cornell.edu/server/api/core/bitstreams/17448f20-41b4-4457-a472-3ff41ab87529/content).

## Discrete charge and energy obstacle to solve first

The present Hall work current is an electric-space Riesz representative:
M j=-C.T H x. Its sampled vector need not have pointwise zero divergence across
element faces. It cannot be treated as a conservative cell electron flux merely
because its weak curl is correct. For a compatible scalar gradient G,
C G=0 implies G.T M j=0, with the appropriate potential/boundary spaces. This
is a dual charge-flux identity; construct and test its actual control-volume
flux interpretation before advancing electron density/thermal energy.

The first implementation now supplies an independent scalar-gradient/current-flux
audit at the axis, conducting wall and periodic seam, followed by a conservative
electron enthalpy transport operator sharing that current and pressure work.
Track which dual fluxes enter primal fluid cells; verify local and boundary
charge conservation, not only a zero global integral. If primal/dual placement
requires a different compatible state space, make the convention explicit and
retain the current Hall baseline as a comparison.

The new current uses K=2pi <r²>g/mu0 and its exact dual face differences. Its
scalar state is nodal with lumped full-domain dual volumes. A physical toroidal
wall trace is required; omitting it inserts a return-current sheet. The isolated
pressure component uses constant wall Btheta/density with zero normal charge
flux, but generally nonzero tangential pressure electric field and local boundary
Poynting flux. It must not silently replace the Hall conducting Ez=0 boundary.

### Next concrete integration gate

Derive an overlap mass map from primal finite-volume ion density to nodal dual
electron density, and a corresponding dual bulk particle flux. In x=r²/2,z,
each primal cell contributes one quarter of its mass to each of its four scalar
vertices, folding periodic endpoints; axis/wall dual volumes account for their
half support. Reconstruct the primal conservative mass flux within each cell
with matching normal face flux and constant divergence (lowest-order H(div)).
Integrate it over dual faces. Prove and independently test that mapping the
primal mass divergence equals the dual flux divergence, including axis/wall/
periodic boundaries. This construction is now implemented and tested locally
against random fluxes and independent affine physical flux integrals.

Then advance electron entropy with bulk dual flux minus J/e while the density
changes, and couple electron pressure/enthalpy and ion energy through an explicit
local conservative total-energy flux. For equal species gammas the total ideal
pressure is fixed by total internal energy, but positive ion/electron partition
still must be verified. Returning a mapped electron-energy change and resetting
cell total energy to fit it is not an acceptable substitute for a derivation.

The ideal moving coupling retains changing n in pe and converges in independent
contact and finite-amplitude compressing simple-wave tests. It keeps total
energy conservative and recovers ion pressure by subtraction. Numerical species
partition and physical-entropy errors remain explicit; this is not a shock
heating closure or entropy inequality guarantee.

The pressure split's fixed-density thermodynamic means alone do not establish
the varying-density chain rule or bulk compression work. Derive those terms,
check exact contact/compression/warm Hall limits and preserve boundary work.
Choose and declare a mutually compatible bounded magnetic/particle condition;
discarding curl-free electric gradients does not turn the physical field into
a conducting-wall field. Keep the eventual true-vacuum/edge model requirement.

### Concrete next warm integration experiment

The current ideal bulk operator is separate from `ElectronPressure`, and its
strong-contact entropy defect is an acceptance case to retain. A candidate
consistent nonideal boundary is fixed Btheta(R)=0, fixed psi wall and free
pressure Ez(R). With zero toroidal wall field, magnetic boundary power is zero
even for varying wall density; this is a mixed magnetic trace, **not** the old
conducting Ez=0 condition. Derive the full Ez electric space/current mass and
verify curl-gradient cancellation, boundary work and independent J1-zero
toroidal Bessel modes before using it. Do not merely remove the density guard
from the existing thermal component and call the full wall consistent.

At fixed n after an ideal half step, the existing pressure/current split can
advance electron U and toroidal B, while the Hall split advances the full
magnetic state. A symmetric composition is a candidate, pending time tests.
Its total-energy increment must come from a **primal face flux** of pressure
Poynting work plus -gamma/(gamma-1) (pe/n) J/e. Reconstruct the divergence-free
dual charge flux within dual cells and integrate enthalpy transport through
primal faces; derive its matching wall traces and consistency independently.
The physical electric gradient must remain in Poynting work even when its
induction curl is zero. Do not substitute the mapped thermal/magnetic energy
change as a local total-energy reset.

Measure both local ion-energy residual and global energy under that split,
with no ion heating expected in the ideal stationary-ion pressure subproblem.
Any local discrete mismatch is a numerical error to converge, not collisional
heat. The independent nonlinear pressure snapshot, pure enthalpy drift and
warm compressible modes should expose it. Only then compose actual moving
warm Hall evolution. Resistivity remains disabled until a justified species
heating model and collision-exchange tests are ready.

## Bounded independent verification sequence

1. Manufactured variable ne, pe and divergence-free J: independently differentiate
   Biermann induction and pressure/enthalpy/Poynting transfers, including axis
   limits. Uniform/barotropic electron pressure supplies exact zero-curl controls.
2. Pure conservative electron thermal drift and compression: verify transport
   with ve, positivity without floors, local pressure work, and exact cancellation
   in total energy. Use fixed-grid time and physical-space refinement, including
   a case with zero bulk velocity but nonzero electron enthalpy transport.
3. Independent homogeneous electron/ion equilibration: analytic relaxation,
   positive temperatures, equal/opposite exchange and correct common temperature.
   Conditional Coulomb coefficients already exist but do not evolve temperature.
4. Warm compressible Hall modes and finite-amplitude thermodynamic response:
   derive the changed continuum oracle and validate species and Maxwell residuals
   before judging the code. Preserve all cold-electron and Hall-free limits.
5. Finite electron inertia: include drift kinetic energy, full convective current
   transport, beta JJ momentum stress and variable-density coefficients from
   the same species reduction. The old fixed-material inertia component is a
   reference limit, not a drop-in nonlinear addition.

Then return to heated-edge compatibility, physical neutral/sheath/kinetic
ordering, true vacuum, spatial coil work, actual gas ends and robust moving
geometry. No conducting atmosphere or density floor may stand in for vacuum.
The original initialized sharp-edge FRC remains a stress case. These gates do
not replace the formation, resolution, topology/lifetime or interactive goals.
