# Fine material density with separate velocity resolution

`MaterialSubcells` is an opt-in one-species block for MaterialSpeciesSet.
Every fine cell retains its actual N, U and physical support. A nested subset
of fine faces carries continuous P1 velocity. No density cut, charge state,
particle, support end or thermal energy is discarded to enlarge the step.
The requested velocity cell count selects nearest cosine-spaced existing faces.

Fine material fractions within a coarse element stay fixed between reactions.
Their physical faces advance with the represented coarse velocity. Geometry is
stored in numpy.longdouble, but this platform supplies only53 mantissa bits:
this currently adds no precision. Crossed or unrepresentable cells reject.

Every block fixes a positive S for its transport substep. Initial input uses
S=max(N_f); source return uses peak species density times a1m reference length.
N=S*N_relative is stored as a factored pair even if an individual physical
amount is below float64 range. Density is evaluated as (N_relative/h)*S, in that
order. No tiny occupied cell is dropped. Internal mass matrix,
force, pressure, Q and thermal state store their physical values divided by S.
This is a coordinate scaling, not a particle threshold or rescaling of physical
N. It prevents a small representable N from underflowing when multiplied by a
kilogram mass. Velocity and temperature are unchanged; source return divides
dimensionless moments by a species scale before multiplying by SI constants.
`physical_thermal` restores physical U for totals. Saved raw states must include
S and explicitly label normalized number and thermal values; physical diagnostic products
below floating-point range can underflow, while the normalized state is retained.

For coarse shapes phi and actual fine density rho=N m/h:

    M_ij = integral(rho phi_i phi_j dx),
    M vdot_j = integral[(p+Q) dphi_j/dx + e Z n E phi_j] dx,
    Udot_f = -(p_f+Q_f) (v_at_right-v_at_left).

The fine velocity jump is the fixed fine/coarse span times the coarse jump.
Analytic polynomial integration builds the scaled tridiagonal consistent mass matrix;
pressure work and forces use identical spans. Electric force uses two-point
Gauss integration over the actual global density union. The returned shape
velocities are also the continuity current. No separate energy repair occurs.

`pressure` retains gamma5/3 and the declared minmod numerical viscosity. Its
compression jump is evaluated on the mechanical mesh; fine sound speed and
mass density enter Q. Numerical heating remains explicit, not physical transport.

`acoustic_timestep` uses cfl/sqrt(lambda_max) from K u=lambda M u, where
K=integral[gamma(p+Q) grad(phi_i) grad(phi_j) dx], plus a deformation bound.
The p part is the linear adiabatic acoustic stiffness; Q is a conservative
positive stress estimate, not a nonlinear stability theorem. MaterialSpeciesSet
also retains its physical electron plasma-frequency restriction. Complete
time/spatial studies must assess both. No implicit stepping or hidden retry.

`restrict_load` maps an exactly gathered fine nodal momentum load to the nested
mechanical space; `fine_velocity` returns represented endpoint velocities for
fine-cell kinetic energy. Projection recovers U=E-K in each actual fine cell;
any negative result rejects. Positive snapshot recovery is not a general theorem.

Only MaterialSpeciesSet supplies the global field and advances this block.
`sample` returns actual fine density/temperature and coarse linear velocity;
the set adds the laboratory frame/current/field. This component does not resolve
collisions, cylindrical fields or walls, and does not establish FRC formation.
