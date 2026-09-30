"""Physical derivatives of regular conservative-moment reconstructions."""
import numpy as np
from scipy.sparse import coo_matrix


def radial_derivative(mesh,points,degree=5,power=0):
    if degree not in (3,5) or power not in (0,1,2) or mesh.nr<degree+1:
        raise ValueError('Supported regular moments and sufficient cells required')
    nodes,w=np.polynomial.legendre.leggauss(degree+2);rows=[];cols=[];data=[]
    for row,point in enumerate(points):
        cell=min(np.searchsorted(mesh.rf,point,side='right')-1,mesh.nr-1)
        for offset in ((degree+1)//2,(degree-1)//2):
            start=min(max(cell-offset,0),mesh.nr-degree-1);ids=np.arange(start,start+degree+1)
            r=mesh.r[ids,None]+.5*mesh.dr*nodes
            scale=.5*(mesh.rf[start+degree+1]**2-mesh.rf[start]**2)
            x=(.5*r*r-.5*point*point)/scale
            weight=.5*mesh.dr*r*w/mesh.radial_measure[ids,None]
            moments=np.stack([np.sum(weight*(r/mesh.radius)**power*x**k,axis=1) for k in range(degree+1)],axis=1)
            target=np.zeros(degree+1)
            target[0]=power*point**(power-1)/mesh.radius**power if power else 0.
            target[1]=(point/mesh.radius)**power*point/scale
            values=.5*np.linalg.solve(moments.T,target)
            rows.extend([row]*(degree+1));cols.extend(ids);data.extend(values)
    return coo_matrix((data,(rows,cols)),shape=(len(points),mesh.nr)).tocsr()


def axial_derivative(mesh,points,degree=5):
    if degree not in (3,5) or mesh.nz<degree+1:raise ValueError('Supported degree and sufficient cells required')
    rows=[];cols=[];data=[]
    for row,point in enumerate(points):
        coordinate=(point-mesh.zf[0])/mesh.dz;cell=min(int(np.floor(coordinate)),mesh.nz-1)
        for offset in ((degree+1)//2,(degree-1)//2):
            ids=np.arange(cell-offset,cell-offset+degree+1);lo=ids-coordinate;hi=lo+1
            moments=np.stack([(hi**(k+1)-lo**(k+1))/(k+1) for k in range(degree+1)],axis=1)
            target=np.zeros(degree+1);target[1]=1/mesh.dz
            values=.5*np.linalg.solve(moments.T,target)
            rows.extend([row]*(degree+1));cols.extend(ids%mesh.nz);data.extend(values)
    return coo_matrix((data,(rows,cols)),shape=(len(points),mesh.nz)).tocsr()
