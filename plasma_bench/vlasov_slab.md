# Conservative electrostatic kinetic slab

## Purpose and equations

An independent planar collisionless component: dt f + v dx f + (qE/m) dv f=0,
dx E=sum(q n)/eps0, E(left)=0. Physical finite ion/electron masses are retained.
This is not the reacting solver or an FRC. A 1D1V distribution is sufficient
for initially factorized Maxwellians and straight Ex: transverse velocities
remain unchanged. Each particle retains perpendicular energy T0*eV, while
parallel pressure and heat evolve directly; no isotropic gamma5/3 closure.

Store g(x,s)=vT*f/n0, s=v/vT, vT=sqrt(T0*eV/m). Uniform finite cells cover the
declared x and s domains. Initial spatial cell averages use the exact slab
overlap; exterior cells are exactly empty. Velocity cell averages integrate
the Gaussian using stable tail differences, without renormalization. A finite
sigma bound omits the explicitly reported Gaussian tail; it is not cutoff-free.

## Components and numerical assumptions

- `advect`: conservative constant-speed transport along each grid line, with
  minmod piecewise-linear reconstruction and a time-centered upwind flux.
  For c>=0, F=c[g+(1-c)*slope/2]; the negative direction mirrors it. Written
  as nonnegative outgoing/remainder contributions, with |c|<=.8. Zero incoming
  boundary flux and measured outgoing flux. No density/pressure floor or repair.
  Minmod causes numerical diffusion and becomes first order at sharp edges;
  there is no physical collision or viscosity term.
- `KineticSpecies`: exact P0 velocity moments through cubic order, including
  cell-width corrections, physical current, directional pressure and central
  parallel heat flux. Zero density alone selects empty velocity outputs.
- `VlasovSlab`: Strang x-half/v-full/x-half transport. Solve Gauss after the
  first half; density is unchanged by acceleration apart from measured velocity
  outflow. Cell-centred Ex drives acceleration. Last Gauss reconstruction supplies
  diagnostics, without charge-neutrality subtraction or field projection.
- `timestep`: explicit spatial/velocity CFL from current field, default .65;
  every actual split stage checks .8. A violated midpoint bound rejects the run,
  not an implicit retry or shortened successful evolution.

## Boundaries and budgets

Zero inflow at remote spatial ends and both velocity limits. Each substep returns
lost N/P/E using P0 velocity integrals for x flux and boundary velocity for v
flux; perpendicular energy travels with number. Boundary losses are not zeroed.
No particle reflection, wall sheath or external circuit is represented.
Gauss E is linear in each x cell; its squared energy is integrated exactly.
Total momentum balance also has electric boundary stress eps0*(ER²-EL²)/2.
The driver records this impulse independently from particle outflow.
`step` also returns the midpoint electric boundary-stress impulse.
Discrete phase transport conserves number with its flux ledger; it does not
enforce exact total energy or momentum. Those defects must decrease with actual
phase-space/time refinement. No post-step renormalization is allowed.

| Dependent | Contract | Important change |
| --- | --- | --- |
| `examples/vlasov_slab.py` | immutable species, tuple of g arrays, explicit losses | phase normalization, moment quadrature, flux or time split |
| saved outputs | physical SI profiles and dimensionless raw moments | zero-current/field constraint would change the physics |

Setting self_consistent=false is the independent free-transport fixture. Its
charge-derived diagnostic E does not act on particles and is excluded from the
conserved free-particle energy. It is not a charged physical experiment.
