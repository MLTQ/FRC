"""Continuous face-flux primitives of mean-preserving cell polynomials."""
import numpy as np
from scipy.sparse import coo_matrix


def difference(edges, periodic=False):
    """Face values to conservative cell divergence in x=r²/2 or z."""
    edges=np.asarray(edges,float);count=len(edges)-1
    rows=np.repeat(np.arange(count),2)
    cols=np.column_stack((np.arange(count),np.arange(1,count+1)))
    if periodic:cols%=count
    data=np.column_stack((-1/np.diff(edges),1/np.diff(edges)))
    return coo_matrix((data.ravel(),(rows,cols.ravel())),
                     shape=(count,count if periodic else count+1)).tocsr()


def primitive_basis(edges, points, degree=5, *, periodic=False):
    """Return cell reconstruction R and flux lift H with dH=R D.

    Each polynomial uses the same two cell-containing stencils as the physical
    conservative density. Its primitive starts at the actual left face value.
    Exact cell moments make the right trace equal the actual right face value.
    """
    edges=np.asarray(edges,float);points=np.asarray(points,float);count=len(edges)-1
    if (edges.ndim!=1 or points.ndim!=1 or degree not in (3,5) or count<degree+1
        or not np.isfinite(edges).all() or not np.isfinite(points).all()
        or np.any(np.diff(edges)<=0) or np.any(points<edges[0]) or np.any(points>edges[-1])):
        raise ValueError('Ordered finite edges and in-domain points with degree 3/5 required')
    if periodic and not np.allclose(np.diff(edges),np.diff(edges)[0],rtol=1e-12,atol=0):
        raise ValueError('Periodic primitive requires uniform axial cells')
    rows=[];cols=[];values=[];integrals=[];left=[]
    for row,point in enumerate(points):
        cell=min(int(np.searchsorted(edges,point,side='right'))-1,count-1)
        left.append(cell)
        for offset in ((degree+1)//2,(degree-1)//2):
            start=cell-offset if periodic else min(max(cell-offset,0),count-degree-1)
            ids=np.arange(start,start+degree+1)
            if periodic:
                lo=edges[0]+ids*(edges[1]-edges[0]);hi=lo+(edges[1]-edges[0])
            else:lo=edges[ids];hi=edges[ids+1]
            origin=.5*(edges[cell]+edges[cell+1]);scale=hi[-1]-lo[0]
            a=(lo-origin)/scale;b=(hi-origin)/scale
            moments=np.stack([(b**(k+1)-a**(k+1))/((k+1)*(b-a))
                              for k in range(degree+1)],axis=1)
            q=(point-origin)/scale;l=(edges[cell]-origin)/scale
            target=np.array([q**k for k in range(degree+1)])
            integral=np.array([scale*(q**(k+1)-l**(k+1))/(k+1)
                               for k in range(degree+1)])
            solved=.5*np.linalg.solve(moments.T,np.column_stack((target,integral)))
            rows.extend([row]*(degree+1));cols.extend(ids%count)
            values.extend(solved[:,0]);integrals.extend(solved[:,1])
    shape=(len(points),count)
    reconstruction=coo_matrix((values,(rows,cols)),shape=shape).tocsr()
    primitive=coo_matrix((integrals,(rows,cols)),shape=shape).tocsr()
    faces=count if periodic else count+1
    anchor=coo_matrix((np.ones(len(points)),(np.arange(len(points)),left)),shape=(len(points),faces))
    return reconstruction,(anchor+primitive@difference(edges,periodic)).tocsr()
