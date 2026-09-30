# Staggered radial state initialization

Independent physical volume quadrature supplies cell means of n, n vz and
n r vtheta. Independent axial quadrature at actual interior radial faces
supplies transverse means of r n vr. Entropies use their declared nodal samples.
This initializes the new DOFs directly from continuum fields; it never labels
an old cell-momentum state as a face state. Initial physical projection error
must be measured independently of subsequent evolution error.
