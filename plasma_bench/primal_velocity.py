"""Mass-weighted regular bulk velocity on the complete dual scalar nodes."""
import numpy as np
from .toroidal_state import radius_squared


def regular_velocity(mapping,conserved):
    mesh=mapping.mesh;rho=conserved[...,0]
    if conserved.shape!=(mesh.nr,mesh.nz,5) or not np.isfinite(conserved).all() or np.any(rho<=0):
        raise ValueError('Finite positive-density conservative moving state required')
    density=mapping.density(rho)
    radial=mapping.density(conserved[...,1]/mesh.radial_centroid[:,None])/density
    rotation=mapping.density(conserved[...,4]/radius_squared(mesh))/density
    axial=mapping.density(conserved[...,2])/density
    # The physical slip-wall condition is explicit; interior values are unchanged.
    radial[-1]=0.
    return np.stack((radial,rotation,axial),axis=-1)
