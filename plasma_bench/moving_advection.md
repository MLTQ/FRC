# Material flux labels and regular axis

Material coordinates are validated before any radial division or stretch use.

`MaterialFlux` preserves each non-axis material node's psi=r²a by storing its
initial a and initial radius. At the axis psi=0 alone cannot determine Bz, so
the regular limit a evolves with the square of the neighboring radial stretch.
Pole nodes use their nearest off-axis angular neighbor; the shared center uses
a sin²-weighted average over the first shell. Uniform radial affine deformation
is exact. Nonuniform axis interpolation is an explicit spatial approximation
that must converge; it is not a continuum frozen-flux theorem.

`coefficients` evaluates current total a; `pullback` is the exact derivative of
this declared discrete rule with respect to radial material coordinates. The
magnetic virtual-work force must include this term in addition to geometric
stiffness derivatives. No field is adjusted to repair energy or equilibrium.

Flux labels at non-axis nodes are invariant under ideal evolution. Continuity
of r²a across shared faces gives a solenoidal weak magnetic field. Preservation
of all interpolated material surfaces/topology is not claimed from nodal
invariance alone. Resistive evolution, remapping and topology change are not
implemented by this ideal-advection component.
