"""Single centered normal-face reconstructions from conservative cell moments."""
import numpy as np
from scipy.sparse import coo_matrix


def axial_faces(mesh,degree=5):
    """Shared periodic face trace centered on an integer face index, not a cell."""
    if degree not in (3,5) or mesh.nz<degree+1:raise ValueError('Degree 3/5 and sufficient axial cells required')
    offsets=np.arange(degree+1)-(degree+1)//2
    moments=np.array([[(float(j+1)**(k+1)-float(j)**(k+1))/(k+1)
                       for k in range(degree+1)] for j in offsets])
    target=np.zeros(degree+1);target[0]=1.;weights=np.linalg.solve(moments.T,target)
    rows=np.repeat(np.arange(mesh.nz),degree+1)
    cols=(np.arange(mesh.nz)[:,None]+offsets)%mesh.nz
    return coo_matrix((np.tile(weights,mesh.nz),(rows,cols.ravel())),shape=(mesh.nz,mesh.nz)).tocsr()


def radial_faces(mesh,degree=5,*,power=1):
    """Centered annular moment stencil, clipped only to the physical interior."""
    if degree not in (3,5) or mesh.nr<degree+1 or power not in (0,1,2):
        raise ValueError('Degree 3/5, regular power 0/1/2 and sufficient radial cells required')
    q,w=np.polynomial.legendre.leggauss(degree+2);rows=[];cols=[];data=[]
    for face,point in enumerate(mesh.rf):
        start=min(max(face-(degree+1)//2,0),mesh.nr-degree-1);ids=np.arange(start,start+degree+1)
        r=mesh.r[ids,None]+.5*mesh.dr*q
        scale=.5*(mesh.rf[start+degree+1]**2-mesh.rf[start]**2);x=(.5*r*r-.5*point*point)/scale
        weight=.5*mesh.dr*r*w/mesh.radial_measure[ids,None]
        moments=np.stack([np.sum(weight*(r/mesh.radius)**power*x**k,axis=1) for k in range(degree+1)],axis=1)
        target=np.zeros(degree+1);target[0]=(point/mesh.radius)**power
        values=np.linalg.solve(moments.T,target)
        rows.extend([face]*(degree+1));cols.extend(ids);data.extend(values)
    return coo_matrix((data,(rows,cols)),shape=(mesh.nr+1,mesh.nr)).tocsr()
