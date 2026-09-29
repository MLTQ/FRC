# Next bounded gate: separate electron and ion thermal transport

The compatible Hall solver now advances nonlinear compressible fields and flow,
but it uses pe=0 and a single represented gas energy. Its prescribed-resistance
heat is a conditional reduced closure, not the original Te=5 eV / Ti=1 eV
argon target. The next step must derive and verify massless electron/ion thermal
transport before extending the moving free surface or claiming formation.

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

Required first implementation: an independent scalar-gradient/current-flux
audit at the axis, conducting wall and periodic seam, followed by a conservative
electron enthalpy transport operator sharing that current and pressure work.
Track which dual fluxes enter primal fluid cells; verify local and boundary
charge conservation, not only a zero global integral. If primal/dual placement
requires a different compatible state space, make the convention explicit and
retain the current Hall baseline as a comparison.

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
