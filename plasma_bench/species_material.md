# Separate material support for each species

`MaterialSpecies` is a planar Euler–Poisson transport component. Each species
has its own material faces x, nodal three-velocity v and cell thermal energy U
(J/m²). Cell particle amounts N (1/m²) and physical masses are fixed in this
reaction-disabled sector. Species fronts evolve with their own endpoint
velocities, not the repeated outer envelope of cell-mean vacuum-fan estimates.
There are no cells or particles outside each support. Empty species are omitted.

For piecewise-linear material velocity, the consistent cell mass matrix is
`m N [[2,1],[1,2]]/6`. Density is N/h. Pressure is `(gamma-1) U/h`.
The semidiscrete equations are

    xdot_j = v_z,j,
    M vdot_z = I.T (p+Q) + integral(e Z n E phi_j dz),
    Udot_i = -(p_i+Q_i) (v_z,i+1-v_z,i).

Incidence I x = x_right-x_left; endpoint forces include no exterior pressure.
Transverse momenta are retained (no transverse force in this planar sector).
Optional compressive numerical Q = rho [C1 c max(-dv,0) + C2 max(-dv,0)^2]
converts its mechanical loss into the same cell thermal energy; it is not
physical viscosity. Defaults C1=C2=0. No energy/temperature reset or floor.

Optional `viscosity_limiter='minmod'` replaces max(-dv,0) in Q by
`max(0, max(-dv_i,0) - h_i min(max(-g_left,0),max(-g_i,0),max(-g_right,0)))`,
where g=dv/h on each species' physical cells. Boundary stencils repeat the
nearest occupied cell, without ghost matter. This removes resolved affine
compression; the residual jump is O(h²) in smooth regions on regular meshes.
Q remains nonnegative and only acts under compression, so its matched work
still heats the actual cell and preserves the semidiscrete energy identity.
Default `none` retains the original explicitly archived stress. This is a
one-dimensional limiter chosen here, inspired by the smooth-compression
criterion of [Caramana et al. (1998)](https://cnls.lanl.gov/~shashkov/papers/ed_vis.pdf),
not a reproduction of that paper's multidimensional viscosity or physical
transport. Its effect on full expansions must be measured before adoption.

The electric field is integrated exactly from the actual piecewise-constant
species charge on the union of all moving cell boundaries. E_left is prescribed;
E_right follows total charge. E is linear on each union interval and constant
through empty intervals out to fixed remote boundaries. Two-point Gauss
integration exactly evaluates E times each species' linear nodal shape.
The represented current is e sum(Z n v); this is the continuity current of
the moving constant-cell-mass map. Its work is the nodal force power, which
cancels the derivative of the integrated interior AND exterior field energy.

Consequently d(K+sum U+field energy)/dt=0 semidiscretely. It is an equation
identity, not an energy repair. RK4 is explicit, acoustic/plasma steps are
disclosed, and actual intermediate/final U and cell widths must stay positive.
Crossed species cells or contact with remote field ends reject without remap.
The consistent mass (not mass lumping) also reproduces affine acceleration for
the quadratic-pressure, uniform-density charged-slab reference.

The staggered pressure/work construction follows the compatible-Lagrangian
principle in [Caramana et al. (1998)](https://cnls.lanl.gov/~shashkov/papers/comp.pdf).
The equations above specify this independent multi-support implementation;
that paper does not establish convergence or stability of this plasma model.

`state`, `step`, `timestep`, `field`, `sample`, and `totals` expose actual
states and budgets to examples.species_material. This component is not adopted
by the formation solver. Atomic reaction transfer between distinct supports,
collisions/heat transport, transverse fields, cylindrical geometry and physical
wall/sheath/outflow physics remain required work. Smooth fluid expansion can
also become kinetically invalid; exact vacuum is not a confinement claim.

Sampling accepts a finite one-dimensional physical coordinate array. All amounts
are per unit transverse area: momentum is kg/(m s), energy J/m² and number 1/m².
