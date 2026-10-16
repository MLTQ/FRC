# Spatial finite-mass species and Maxwell fields

`AtomicMaxwell` extends the charge-continuous longitudinal model with both
transverse electric/magnetic components and all three components of the actual
Lorentz force. It uses the full species number/momentum/energy state and physical
electron mass, epsilon0 and light speed. No reduced light speed or quasineutral
projection is used. The fixed laboratory interval has reflecting particle walls
and perfectly conducting electromagnetic ends (Et=0). These are verification
boundaries, not a quartz tube, sheath, open end or tracked plasma interface.

## Spatial equations

Coordinate z follows the inherited axial transport. Transverse E=(Ex,Ey) is
stored on faces; B=(Bx,By) at cell centres. Bz is a constant input, so div B=0
identically in 1D. Longitudinal E follows the actual charge through the inherited
Gauss solve and numerical particle-flux current; see atomic_electrostatic.md.
The transverse equations are

    dEx/dt = -c² dBy/dz - Jx/epsilon0
    dEy/dt = +c² dBx/dz - Jy/epsilon0
    dBx/dt = +dEy/dz; dBy/dt = -dEx/dz.

Cell transverse J=e Nref sum Z_s f_s v_s comes from actual extensive momenta.
Arithmetic averaging scatters it to interior E faces. The transpose averaging
gathers E to cell forces and physical work q n v·E. This makes source work cancel
field work in the semidiscrete face/cell energy. The centred curl operators are
opposite adjoints under the same quadrature and zero tangential boundary field;
Poynting flux through these ends vanishes. No field-energy residual is fed back.

Species momentum gains q n(E_transverse+v cross B) in addition to the inherited
longitudinal electric force. Magnetic force does no laboratory work. Stored
fixed-frame energy also receives minus U dot force. Both empty matter and empty
species remain exactly zero. The equations are the nonrelativistic five-moment
multifluid–Maxwell model described by [Wang et al.](https://arxiv.org/abs/1909.04125);
the staggered spatial discretization and longitudinal numerical-flux work above
are specified here independently, rather than claimed identical to that paper.

## Integration and measured limits

RK4 advances the entire spatial/field force system together; no frozen-source
splitting error is hidden inside this operator. Its explicit step resolves the
actual acoustic CFL, c dt/dz, local plasma frequency and gyrofrequency. RK4 is
chosen over SSP-RK2 because a centred undamped Maxwell curl has imaginary modes;
RK2 has no stable imaginary-axis interval. RK4 still has finite phase/energy
error. Local reaction/collision halves are optional and retain their own accuracy
limits; the wave runs explicitly disable them.

`totals` measures particle+binding+photon-input+longitudinal+transverse field
energy. Constant Bz energy is omitted as time independent. It also reports the
energy using exact piecewise-linear E integration, to expose the face-lumped
quadrature difference. `last_diagnostic` measures semidiscrete work cancellation;
it never repairs any state. Wall stress can exchange momentum with the domain;
no closed-system momentum claim is made. A discrete full momentum/stress balance
and spatial reaction/collision coupling still need measurement.

HLLC/HLL/minmod particle diffusion is inherited; centred Maxwell curls add no
spatial damping. The Eulerian vacuum-tail limitation remains, so this component
must not replace the formation model's tracked exact exterior. Axisymmetric
geometry, real coil/open-end boundaries, physically valid transport and a moving
interface remain required. Wave dispersion is not formation or confinement.
