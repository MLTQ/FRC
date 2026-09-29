"""Smooth physical velocity from actual conservative annular momentum averages."""
import numpy as np
from scipy.sparse import coo_matrix


def radial_moments(mesh,points,degree=5,*,power=0):
    if degree not in (3,5) or mesh.nr<degree+1 or power not in (0,1,2):raise ValueError('Degree 3/5, regular power 0/1/2 and sufficient cells required')
    nodes,weights=np.polynomial.legendre.leggauss(degree+2)
    rows=[];cols=[];data=[]
    for row,point in enumerate(points):
        cell=min(np.searchsorted(mesh.rf,point,side='right')-1,mesh.nr-1)
        for offset in ((degree+1)//2,(degree-1)//2):
            start=min(max(cell-offset,0),mesh.nr-degree-1);ids=np.arange(start,start+degree+1)
            r=mesh.r[ids,None]+.5*mesh.dr*nodes
            scale=.5*(mesh.rf[start+degree+1]**2-mesh.rf[start]**2)
            x=(.5*r*r-.5*point*point)/scale
            weight=.5*mesh.dr*r*weights/mesh.radial_measure[ids,None]
            moments=np.stack([np.sum(weight*(r/mesh.radius)**power*x**k,axis=1) for k in range(degree+1)],axis=1)
            target=np.zeros(degree+1);target[0]=(point/mesh.radius)**power
            values=.5*np.linalg.solve(moments.T,target)
            rows.extend([row]*(degree+1));cols.extend(ids);data.extend(values)
    return coo_matrix((data,(rows,cols)),shape=(len(points),mesh.nr)).tocsr()


def axial_moments(mesh,points,degree=5):
    if degree not in (3,5) or mesh.nz<degree+1:raise ValueError('Degree 3/5 and sufficient periodic cells required')
    rows=[];cols=[];data=[]
    for row,point in enumerate(points):
        coordinate=(point-mesh.zf[0])/mesh.dz;cell=min(int(np.floor(coordinate)),mesh.nz-1)
        for offset in ((degree+1)//2,(degree-1)//2):
            ids=np.arange(cell-offset,cell-offset+degree+1);lo=ids-coordinate;hi=lo+1
            moments=np.stack([(hi**(k+1)-lo**(k+1))/(k+1) for k in range(degree+1)],axis=1)
            target=np.zeros(degree+1);target[0]=1.
            values=.5*np.linalg.solve(moments.T,target)
            rows.extend([row]*(degree+1));cols.extend(ids%mesh.nz);data.extend(values)
    return coo_matrix((data,(rows,cols)),shape=(len(points),mesh.nz)).tocsr()


def velocity(geometry,conserved,*,degree=5):
    """Sample only; does not replace continuity, density state or fluid energy."""
    mesh=geometry.mesh;u=np.asarray(conserved)
    if u.shape!=(mesh.nr,mesh.nz,5) or not np.isfinite(u).all() or np.any(u[...,0]<=0):raise ValueError('Finite positive-density conservative state required')
    axial=axial_moments(mesh,geometry.axial_points,degree)
    radial={power:radial_moments(mesh,geometry.radial_points,degree,power=power) for power in (0,1,2)}
    sampled=[geometry.tensor(radial[power],u[...,column],axial) for column,power in ((0,0),(1,1),(2,0),(4,2))]
    rho,mr,mz,angular=sampled
    if np.any(rho<=0) or not np.isfinite(sampled).all():raise ValueError('Nonpositive reconstructed density; no floor or fallback')
    return np.stack((mr/rho,angular/(rho*geometry.r),mz/rho),axis=-1)
