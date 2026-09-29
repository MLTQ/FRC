"""Shared-flux smooth transport of dual conservative entropy and density."""
from functools import lru_cache
import numpy as np
from .annular_polynomial import axial_samples
from .mhd_state import NonphysicalState


@lru_cache(maxsize=32)
def radial_faces(edges):
    edges=np.asarray(edges);count=len(edges)-1;matrix=np.zeros((count+1,count))
    matrix[0,0]=matrix[-1,-1]=1.
    for face in range(1,count):
        start=min(max(face-2,0),count-4);ids=np.arange(start,start+4)
        scale=edges[start+4]-edges[start];lo=(edges[ids]-edges[face])/scale;hi=(edges[ids+1]-edges[face])/scale
        moments=np.stack([(hi**(p+1)-lo**(p+1))/((p+1)*(hi-lo)) for p in range(4)],axis=1)
        matrix[face,ids]=np.linalg.solve(moments.T,np.array([1.,0.,0.,0.]))
    matrix.setflags(write=False);return matrix


def entropy_rate(mapping,density,entropy,radial,axial):
    geo=mapping.geometry;n=geo.scalar(density,'number density');u=geo.scalar(entropy,'entropy proxy')
    radial=np.asarray(radial);axial=np.asarray(axial)
    if radial.shape!=(mapping.mesh.nr+2,mapping.mesh.nz) or axial.shape!=u.shape or not np.isfinite(radial).all() or not np.isfinite(axial).all() or np.any(radial[[0,-1]]!=0):
        raise ValueError('Finite shared dual particle flux with impermeable radial boundary required')
    matrix=radial_faces(tuple(geo.dual_x_faces));nr,ur=matrix@n,matrix@u
    nz,uz=axial_samples(n,1.),axial_samples(u,1.)
    if any(not np.isfinite(v).all() or np.any(v<=0) for v in (n,u,nr,ur,nz,uz)):
        raise NonphysicalState('Nonpositive smooth reconstructed density/entropy; no repair or fallback')
    return -geo.divergence(radial*(ur/nr),axial*(uz/nz))
