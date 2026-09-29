# Electron pressure and entropy split component

At fixed positive n, this component pairs pressure-driven toroidal induction
with adiabatic entropy transport by ve=-J/(n e). It does not include bulk fluid
motion, Hall J cross B, resistivity, collisions, conduction or electron inertia.
It is not connected to HallMHD's conducting-wall fluid state yet.

Let a=(n/nref)^(gamma-1), w=p/(n a), U=n w, A=a/(gamma-1). Electron internal
energy is A U. The continuum equations for this split are
Udot=div(w J/e) and Epressure=-grad(p/n)/e-w grad(A)/e=-grad(p)/(n e).
The same thermodynamic edge w is used in electron transport and electric field:
the reciprocal logarithmic mean of the reciprocal nodal w values.
The current's exact dual divergence cancels the pure gradient work. Thus
magnetic and electron thermal energy exchange exactly in the assembled pairing;
constant w gives exactly curl-free (adiabatic barotropic) electric pressure.
For constant temperature, w=T/a and this mean makes w_face delta(A) exactly
T delta(log n), so isothermal barotropic induction is also curl-free.
The identity w_face delta(1/w)=-delta(log w), together with local charge
conservation, preserves physical logarithmic entropy semidiscretely.
For induction, the pure gradient is omitted using C grad=0. Within each closed
A circulation a local constant w is subtracted, another exact gradient null
identity. This avoids catastrophic axis cancellation of large pressure-gradient
terms; the full electric field is still returned for physical work/flux checks.
It is an algebraic evaluation of the same curl, not an axis correction or floor.
Fixed density and derived thermodynamic coefficients are immutable arrays.

The full-domain dual nodes carry scalar values and lumped volumes. This is an
explicit change in placement for the isolated thermal component, not automatic
conversion of existing primal ion/fluid energy. Density and wall Btheta are
fixed; Btheta must be constant in z. For nonzero Btheta, wall density must also
be constant. For **zero Btheta**, arbitrary positive wall density is allowed:
normal wall current and pointwise magnetic boundary power are both zero.
Ez at the wall is a pressure-gradient field, generally nonzero. For constant
nonzero Btheta/density its integrated Poynting input still vanishes. This trace
boundary is **not** HallMHD's conducting Ez=0 boundary. General particle inflow
and sheath/reservoir boundary thermodynamics remain unimplemented.

The symmetric entropy-conservative transport has no monotonicity theorem.
Positive states are required and failure raises without flooring. Physical
logarithmic entropy is diagnosed independently of conserved linear proxy U.
The temporal logarithmic mean in electron_midpoint preserves its finite-step
identity too; local thermodynamic/field errors still must converge. Energy is not
proof of thermodynamic accuracy. The arbitrary positive density scale changes
the proxy units, not physical p, induction, thermal work or particle transport.
