"""Cell-conservative transverse cubic polynomials for face quadrature."""
import numpy as np


def radial_cell_matrix(mesh,fraction,*,power=0):
    """Evaluate each cell's mean-preserving polynomial at its radial fraction."""
    if power not in (0,1,2) or not np.isfinite(fraction) or not 0<=fraction<=1:
        raise ValueError('Regular power 0/1/2 and a fraction within the cell required')
    matrix=np.zeros((mesh.nr,mesh.nr));nodes,weights=np.polynomial.legendre.leggauss(5)
    for cell in range(mesh.nr):
        r=mesh.rf[cell]+fraction*mesh.dr
        # Average the two mirrored cubic stencils; clip to physical interior
        # cells at the axis/wall. Both contain this cell and recover its mean.
        for start in (min(max(cell-2,0),mesh.nr-4),min(max(cell-1,0),mesh.nr-4)):
            ids=np.arange(start,start+4);scale=.5*(mesh.rf[start+4]**2-mesh.rf[start]**2)
            rr=mesh.r[ids,None]+.5*mesh.dr*nodes[None]
            x=(.5*rr*rr-.5*r*r)/scale
            volume=.5*mesh.dr*rr*weights[None]/mesh.radial_measure[ids,None]
            basis=(rr/mesh.radius)**power
            moments=np.stack([np.sum(volume*basis*x**p,axis=1) for p in range(4)],axis=1)
            target=np.array([(r/mesh.radius)**power,0.,0.,0.])
            matrix[cell,ids]+=.5*np.linalg.solve(moments.T,target)
    return matrix


def axial_cell_samples(value,fraction):
    """Centered mean of mirrored, mean-preserving cubic cell polynomials."""
    if not np.isfinite(fraction) or not 0<=fraction<=1:raise ValueError('Fraction inside cell required')
    result=np.zeros_like(value,dtype=float)
    for offsets in (np.arange(-2,2),np.arange(-1,3)):
        moments=np.array([[(float(k+1)-fraction)**(p+1)/(p+1)-(float(k)-fraction)**(p+1)/(p+1)
                           for p in range(4)] for k in offsets])
        weights=np.linalg.solve(moments.T,np.array([1.,0.,0.,0.]))
        result+=.5*sum(w*np.roll(value,-int(k),axis=1) for w,k in zip(weights,offsets))
    return result
