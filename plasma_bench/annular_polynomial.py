"""Cubic interpolation of true annular cell moments for smooth reference fluxes."""
import numpy as np


def radial_matrix(mesh,radii,*,power=0):
    """Map cell averages of r^power P3(r²/2) to supplied physical radii."""
    radii=np.asarray(radii,float)
    if radii.ndim!=1 or not np.isfinite(radii).all() or np.any(radii<0) or np.any(radii>mesh.radius) or power not in (0,1,2):
        raise ValueError('Radii inside the cylinder and regular radial power 0, 1 or 2 required')
    centers=.25*(mesh.rf[:-1]**2+mesh.rf[1:]**2);matrix=np.zeros((len(radii),mesh.nr))
    nodes,weights=np.polynomial.legendre.leggauss(5)
    for row,r in enumerate(radii):
        start=min(max(int(np.searchsorted(centers,.5*r*r))-2,0),mesh.nr-4);ids=np.arange(start,start+4)
        xscale=.5*(mesh.rf[start+4]**2-mesh.rf[start]**2)
        rr=mesh.r[ids,None]+.5*mesh.dr*nodes[None]
        x=(.5*rr*rr-.5*r*r)/xscale
        volume=.5*mesh.dr*rr*weights[None]/mesh.radial_measure[ids,None]
        basis=(rr/mesh.radius)**power
        moments=np.stack([np.sum(volume*basis*x**p,axis=1) for p in range(4)],axis=1)
        target=np.array([(r/mesh.radius)**power,0.,0.,0.])
        matrix[row,ids]=np.linalg.solve(moments.T,target)
    return matrix


def axial_samples(value,fraction):
    """Cubic reconstruction from periodic axial cell averages at zface_j+fraction dz."""
    if not np.isfinite(fraction) or not 0<=fraction<=1:raise ValueError('Axial fraction within a cell required')
    offsets=np.arange(-2,2) if fraction<.5 else np.arange(-1,3)
    moments=np.array([[(float(k+1)-fraction)**(p+1)/(p+1)-(float(k)-fraction)**(p+1)/(p+1)
                       for p in range(4)] for k in offsets])
    weights=np.linalg.solve(moments.T,np.array([1.,0.,0.,0.]))
    return sum(w*np.roll(value,-int(k),axis=1) for w,k in zip(weights,offsets))


def normal_momentum_flux(mesh,conserved):
    """Candidate high-order continuity flux; explicit zero radial particle boundary."""
    if conserved.shape!=(mesh.nr,mesh.nz,5) or not np.isfinite(conserved).all():
        raise ValueError('Finite conservative moving state required')
    radial=radial_matrix(mesh,mesh.rf,power=1)@conserved[...,1]
    radial[[0,-1]]=0.
    axial=axial_samples(conserved[...,2],0.)
    axial=np.concatenate((axial,axial[:,:1]),axis=1)
    return 2*np.pi*mesh.rf[:,None]*mesh.dz*radial,2*np.pi*mesh.radial_measure[:,None]*axial
