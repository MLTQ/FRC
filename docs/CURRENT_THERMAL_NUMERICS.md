# Finite-mass current and species thermal reference

`CurrentThermal` is an experimental, bounded, **fixed-density and constrained
zero-barycentric-flow** subsystem. It includes physical electron mass and full
nonlinear current inertia. It does not yet integrate moving finite-mass plasma,
vacuum, gas ends, prescribed coils, collisions or kinetic pressure closure.
`WarmHallMHD` and guarded Stage2 have not been silently switched to this model.

## Parent species and the explicit constraint

For quasineutral singly charged species define

    M=mi+me, rho=M n, alpha=me/M, beta=mi/M,
    delta=beta-alpha, mu=mi me/M, w=J/(n e),
    vi=v+alpha w, ve=v-beta w, J=curl B/mu0.

Here me is the physical electron mass; `model='hall'` consistently sets it to
zero as a distinct comparison. There is no adjustable effective electron mass.
With isotropic adiabatic pressures and no collisions, parent momentum gives

    rho (partial_t+v dot grad) v
       = J cross B - grad(pi+pe) - div(mu n w w)

    E+v cross B = delta w cross B + (alpha grad pi-beta grad pe)/(n e)
       + mu/e [partial_t w + (v dot grad)w + (w dot grad)v
               - delta(w dot grad)w].

Both species obey `partial_t n + div(n vs)=0` because div J=0, and
`partial_t ps + vs dot grad ps + gamma_s ps div vs=0`.
Displacement current, charge separation, pressure anisotropy, nongyrotropic
stress and kinetic transport are absent. Numerical verification of this
reduction cannot establish that those omissions are physical for an FRC edge.

The implemented subproblem sets v=0 and holds positive n(x) fixed. This requires
an external bulk constraint force

    R = mu n (w dot grad)w - J cross B + grad(pi+pe),

distributed as beta R on ions and alpha R on electrons. It cancels in the
difference equation and does zero total work because beta vi+alpha ve=0.
That constraint is exposed, not claimed to be a physical unforced equilibrium.

The species framework follows the [NRL formulary, pp.36–39](https://www.govinfo.gov/content/pkg/GOVPUB-D210-PURL-LPS118567/pdf/GOVPUB-D210-PURL-LPS118567.pdf).
Finite-mass/current-energy context is also discussed by
[Martin, sections 2.2–2.3](https://ecommons.cornell.edu/server/api/core/bitstreams/17448f20-41b4-4457-a472-3ff41ab87529/content)
and [Andrés et al.](https://arxiv.org/pdf/1410.6518).
The exact mass fractions, constraint and numerical construction here are this
project's explicit derivation; those sources do not validate this code or its
experimental regime, and their additional reductions are not automatically used.

## Current inertia and full energy flux

For this fixed-n subproblem let

    B* = B + curl(mu J/(n e²)) = B + mu/e curl w
    KJ = mu J²/(2 n e²) = mu n |w|²/2
    epsilon_s=ps/(gamma_s-1), Hs=epsilon_s+ps.

Then the equivalent Ohm equation is

    E = delta w cross B* + mu/(n e²) partial_t J
        - mu delta/(2e) grad |w|²
        + (alpha grad pi-beta grad pe)/(n e).

Retaining only partial_t J would omit nonlinear convective inertia. The
implementation retains B* **and** the scalar kinetic gradient. Its energy law is

    partial_t [B²/(2mu0)+KJ+epsilon_e+epsilon_i] + div F = 0
    F = E cross B/mu0 + w [alpha Hi-beta He-delta KJ].

The kinetic scalar enters both full E and the drift energy flux. A zero curl
does not make its Poynting work dispensable. Species entropies
`Us=ps/(n/nref)^(gamma_s-1)` obey

    partial_t Ue = beta div[(Ue/n) J/e]
    partial_t Ui = -alpha div[(Ui/n) J/e].

## Discrete update and boundaries

The periodic filled cylinder retains fixed axis/wall psi, Btheta(R)=0 and free
Ez(R). It is a bounded magnetic trace, not a conducting electric wall or vacuum.
Compatible magnetic/electric spaces in x=r²/2,z give C=-curl and weak current
`j=-M_electric^-1 C.T H b`. Positive Q1 density is sampled on common quadrature;
`M_inertia` represents integral mu |J|²/(n e²). Density and inertia coefficients
are immutable. These matrices describe actual physical current energy.

One nonlinear solve uses magnetic/current midpoint and temporal logarithmic
means of each Us. Spatial entropy means retain pressure-work pairing. Hall
transport is skew with delta B*. The mixed solve uses current increments:

    M_electric dj_mid + dt/2 C.T H C E_reduced = 0
    M_electric E_reduced - 2/dt M_inertia dj_mid
       = X(delta B*) (jold+dj_mid) + M_electric Ep_reduced.

Pure gradients are removed from the induction solve and restored for full E.
Pressure includes beta electron minus alpha ion pressure electric fields.
Sparse linear and nonlinear residuals are checked explicitly. The current
change is evaluated directly as `j(db/2)`, not by subtracting two background
currents. A separate Ampere/Faraday backward residual checks the stored state
on the absolute matrix-product/state scale; it cannot replace the relative
small-signal Ohm residual. This distinction resolves stationary-bias cancellation
without changing tolerances or suppressing any physical current.
Nonpositive species, nonfinite fields or unconverged trials reject atomically.
There is no hidden retry, temperature floor, effective-mass adjustment or heat.

Face-integrated Poynting, species enthalpy and drift-energy flux are returned
in watts. The kinetic potential/flux uses a stated nodal current interpolation,
including normal-current and azimuthal boundary traces. Its **local** accuracy
is measured; a global work identity is not an exact cell energy theorem.
Verification evolves a separate conservative total-cell energy solely from
these face powers and compares it with magnetic/current/species component
energies. It also rejects nonpositive conservative ion partition. No reset
enforces agreement between these two energy accounts.

## Exact finite-amplitude wave

Take q a J1 zero divided by R, k=2pi m/L, K=sqrt(q²+k²), h=+/-1,
and a=A J1(qr) exp[i(kz-omega t)]. With uniform n and pressures,

    b=(-i k a, h K a, A q J0(qr) exp[i(kz-omega t)])
    B=B0 ez+Re b, curl b=h K b
    omega = delta h B0 k K / [mu0 n e (1+mu K²/(mu0 n e²))].

J is parallel to the perturbation B*, so the nonlinear cross term vanishes at
finite amplitude. The independent reference differentiates cylindrical drift
velocity analytically, exposes the constraint reaction, and computes full E
from parent species acceleration. Tests substitute it into both species
momentum equations, zero constraint work, species kinetic energy, Faraday and
the **local** conservative total-energy equation. Both helicities are retained.
Actual field/current and full face-power errors are then measured, not inferred
from global conservation. CURRENT_THERMAL_VERIFICATION records measured limits.

## Required moving extension

For general v,n the total flux is

    F = (Kbulk+KJ+Hi+He) v + E cross B/mu0
        + w [mu n v dot w - delta KJ + alpha Hi-beta He],
    Kbulk=rho |v|²/2.

The bulk momentum drift stress, changing 1/n current energy, bulk/current
advection, velocity-gradient terms and both species bulk entropy fluxes are
required together. Updating only the current mass matrix as density moves is
not sufficient. Derive a conservative coupled stage before inserting this
fixed component in the existing ideal/Hall splitting. Verify the full finite-
mass compressible wave spectrum, variable-density force/energy transport and
local pressure partition independently, retaining maximum/axis errors and
failed trials. Collision heating/exchange, conduction/ionization, the original
heated free edge, true vacuum, real ends/coils and physical validity remain
separate requirements. None is discharged by this constrained reference.

The subsequent MOVING_CURRENT_NUMERICS/VERIFICATION checkpoint derives and
independently checks the exact moving continuum contract and instantaneous
induction rate, including n_t and curl v. It leaves the actual conservative
moving fluid time update as the next integration gate; this fixed-n module
and its frozen measurements remain unchanged.
