# Instantaneous interface-energy diagnostic

The local limiter study separates improved central shock profiles from an
unresolved full-space energy error near the transparent gas wall. The following
identity provides a direct diagnostic. Its independent perturbation check and
measured values are in the [local-limiter study](LOCAL_LIMITER_VERIFICATION.md);
it is not a proposed thermal-energy correction.

At a fixed accepted state, write the fluid total-energy outward rate as Pout,
the full induced unknown vector as u, and the compatible cell magnetic field
as Bc. Vacuum stationarity gives (K u)v=0. If g is the total corner-flux rate
on active plasma nodes, then induced active rate is g minus prescribed-source
rate. Differentiate the independently defined response energy and subtract
the source power used by the solver:

```
Eresponse = sum(Ecell V) - sum(Bc,total² V / (2 mu0)) + pi u^T K u / mu0
Psource = -(2 pi / mu0) sum_active[(K u) dpsi_external/dt]

Rdot = dEresponse/dt - Psource
     = -Pout
       - sum_cells[Bc,total · cell_field(g) V / mu0]
       + (2 pi / mu0) sum_active[(K u) g]
```

The source terms cancel algebraically; all-space induced energy still includes
vacuum/multipole contributions through K. The vacuum rows contribute no work,
but cannot be removed from the magnetic completion or from K. The corner-rate
array must be embedded in the field grid with the existing active-node mapping.
This identity assumes the same fixed fluid/vacuum partition as the reference
solver; a moving interface will require its own derivation.

Since g=-r Etheta and (K u)=mu0 Jtheta Adual on active nodes, the last term is
-2 pi sum(Jtheta Etheta r Adual). It includes ideal magnetic work, physical
resistivity and declared interface numerical diffusion. It is not only Joule
heat. It must not be added to the total fluid energy a second time.

For selected local versus original high-order fluxes at the **same state**:

```
delta Rdot = -delta Pout
             - sum[Bc,total · cell_field(delta g) V / mu0]
             + (2 pi / mu0) sum_active[(K u) delta g].
```

Only changed boundary fluxes and corner stencils contribute to this difference.
It can isolate the direct limiter/interface mismatch without comparing two
already diverged trajectories or attributing a global energy cancellation to
interior accuracy. The high, low and locally selected rates should all be
recorded at common states, alongside their masks and the cumulative measured
response residual.

For an independent check, perturb the state in the selected conservative and
corner-flux direction, complete its induced vacuum, and differentiate the
independently integrated response energy. For the constant-source shock, this
energy is quadratic along the perturbation; a symmetric difference cancels the
quadratic term. Vary perturbation size to distinguish roundoff from a mapping or
factor-of-two error. Positive-state perturbations must be checked explicitly.

Next numerical work should use this identity to distinguish spatial interface
flux mismatch from RK truncation and from indirect changes in the evolving wall
layer. Any proposed correction must change consistent conservative/CT fluxes
and pass independent boundary-work and refinement checks. Do not reset thermal
energy, hide the residual with a global normalization, or interpret a small
vacuum linear-solve residual as a physical energy-accuracy certificate.
