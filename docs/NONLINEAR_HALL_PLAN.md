# Next integration: nonlinear Hall work on the compressible toroidal state

The bounded compressible Hall milestone is now implemented and independently
verified in COMPATIBLE_HALL_NUMERICS and COMPATIBLE_HALL_VERIFICATION. This
document retains its integration requirements and historical failed-axis gate.
TWO_TEMPERATURE_PLAN is the next bounded task; full formation scope is unchanged.

## Historical integration gate and completed correction

The new finite-pressure oracle and nonlinear variable-density manufactured
reference are implemented (COMPRESSIBLE_HALL_REFERENCE). An experimental
common-quadrature Hall operator now cancels assembled magnetic work and has
explicit Poynting transport, but **fails the regular-axis consistency gate**:
the first-annulus gdot error approaches 21.46% through nr=128. It is not enabled
in ToroidalMHD. The 42 completed instantaneous cases and source/archive evidence
are in HALL_WORK_VERIFICATION; equations are in HALL_WORK_NUMERICS.

The immediate next action is to reconcile regular polynomial moments, current
and electric sampling with the magnetic energy convention. Use exact local
polynomials and the retained full nonlinear snapshot, not a fitted correction
for this one case. Any energy-quadrature change must be explicit in state
pressure recovery, fluxes and diagnostics. Preserve the existing verified MHD
path and rerun relevant equilibrium/transport tests. Do not enable Hall merely
because its global work cancels or a shrinking-axis-volume norm converges.

The compatible successor now passes that axis gate with an explicit new energy
convention and a physical commuting curl. Local transfer, nonlinear midpoint,
actual compressible waves and variable-density evolution are verified in 63
completed cases; the old failed operator is preserved. First-annulus error is
8.30089e-5 at nr128, and fixed-grid wave/nonlinear time errors fall near fourfold
per halving. Physical magnetic representation remains first order. This satisfies
the bounded cold-electron integration milestone, not thermal or vacuum closure.

The original integrator requirement was to verify local energy transfer and a
stable time integrator. Purely adding
nondissipative Hall to the existing SSPRK2 step is not stable for imaginary
wave eigenvalues. Actual nonlinear wave propagation, mesh/time refinement,
variable-density transfer and finite-amplitude controls are now implemented;
see the newer report for their measured errors and durations.

## Current and electromagnetic work

For singly charged massless electrons with prescribed fully ionized species,

    ne=rho/mi, E_H=(J cross B)/(ne e)
    E_H,r=(Jtheta Bz-Jz Btheta)/(ne e)
    E_H,theta=(Jz Br-Jr Bz)/(ne e)
    E_H,z=(Jr Btheta-Jtheta Br)/(ne e).

These terms enter both poloidal CT and toroidal Faraday, while their Poynting
flux E_H cross B/mu0 enters total energy. They do not supply a Joule heat source.
The fluid magnetic force already uses full Maxwell stress; it must not receive
a second Hall force. A common compatible current/field quadrature must establish
J dot E_H=0 in the assembled work pairing, with explicit boundary terms and
numerical dissipation. Merely replacing fluid velocity by v-J/(ne e) in one
induction equation omits required terms and is insufficient.

Positive carrier density inside this bounded reference is an assumption, not
a vacuum regularization. The free vacuum region still has no mass/current. A
later interface cannot divide by zero density or use a floor to mimic carriers.
Electron pressure adds -(grad pe)/(ne e); it is not generally curl-free for
variable ne. A nonlinear massless-Hall checkpoint may use an explicitly stated
barotropic electron-pressure reduction, but cannot claim the full thermal goal.
The implemented reference chooses explicitly cold electrons pe=0 and adiabatic
ion pressure; a warm barotropic extension also needs electron enthalpy/work
accounting even when its induction curl vanishes.

## Independent compressible cylindrical wave reference

The earlier cylindrical reference is incompressible. At finite sound speed,
use a separate continuum reduction with radial J1(qr) transverse variables and
J0(qr) scalar/axial variables; qR is a J1 zero and k=2pi ell/L. Let
Br=-ik a J1, Bz=q a J0, Btheta=t J1, vr=u J1, vtheta=w J1, vz=v J0,
and pressure perturbation p1=p J0. Linearizing massless Hall about uniform
rho, p0, B0 ez gives six amplitude equations:

    rho u_dot = q p + B0 (q²+k²) a/mu0
    rho w_dot = ik B0 t/mu0
    rho v_dot = -ik p
    p_dot = -gamma p0 (q u+ik v)
    a_dot = -B0 u - ik B0 t/(mu0 ne e)
    t_dot = ik B0 w - ik B0 (q²+k²) a/(mu0 ne e).

These follow from continuum momentum, pressure and Faraday and contain both
compressible and torsional branches. Independently check their species/Maxwell
residuals and limiting dispersion before using them as a numerical oracle.
For eta=0 their homogeneous traces satisfy the bounded conducting wall.
Finite eta needs its own compatible current boundary: the earlier J1 toroidal
trace has nonzero Jz at that wall and must not be mislabeled resistively
conducting. Use independent periodic/analytic references for damping if needed.

## Nonlinear and finite-mass gates

1. Check both helicities/propagation signs, sound-speed and Hall-length limits,
   actual mesh/time refinement, global and local energy transfer, and divergence.
2. Exercise finite-amplitude Hall redistribution with independent manufactured
   residuals or a separately converged reference, retaining the original
   nonlinear rotation/pinch/contact tests and zero-Hall reduction.
3. Add electron/ion thermal equations, collision exchange and pressure work.
   Finite electron mass requires drift kinetic energy, beta J J momentum stress,
   convective current transport and variable-density coefficients from the same
   parent equations; the fixed-current component cannot simply be bolted on.
4. Rejoin the exact-vacuum/moving-boundary path with spatial coil source work,
   physical edge/neutral/kinetic ordering, real ends and converged formation.

The original low-pressure argon target, heated-edge incompatibility and interface
requirements remain in EXTENDED_PLASMA_PLAN. This bounded Hall work is the next
integration milestone, not a replacement completion criterion.
