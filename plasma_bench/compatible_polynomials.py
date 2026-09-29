"""Local nodal polynomials and interval-average edge bases with commuting d."""
import numpy as np
from scipy.sparse import coo_matrix


def bases(edges,points,degree=5,*,periodic=False):
    """Return nodal N, dN/dx, edge-average E, and nodal-difference D."""
    edges=np.asarray(edges,float);points=np.asarray(points,float)
    count=len(edges)-1
    if degree not in (1,3,5) or count<degree+int(periodic):
        raise ValueError('Degree 1/3/5 and enough independent grid nodes required')
    if edges.ndim!=1 or np.any(np.diff(edges)<=0) or not np.isfinite(edges).all() or points.ndim!=1 or not np.isfinite(points).all():
        raise ValueError('Ordered finite edges and finite one-dimensional points required')
    if np.any(points<edges[0]) or np.any(points>edges[-1]):raise ValueError('Points must lie inside the domain')
    if periodic and not np.allclose(np.diff(edges),np.diff(edges)[0],rtol=1e-13,atol=0.):
        raise ValueError('Periodic basis currently requires uniform intervals')
    n=count if periodic else count+1
    nr=[];nc=[];nv=[];dv=[];er=[];ec=[];ev=[]
    for row,point in enumerate(points):
        cell=min(np.searchsorted(edges,point,side='right')-1,count-1)
        start=cell-(degree-1)//2
        if not periodic:start=min(max(start,0),count-degree)
        ids=np.arange(start,start+degree+1)
        coords=edges[0]+ids*(edges[1]-edges[0]) if periodic else edges[ids]
        scale=coords[-1]-coords[0];t=(coords-point)/scale
        vandermonde=np.vander(t,N=degree+1,increasing=True).T
        target=np.zeros(degree+1);target[0]=1.
        value=np.linalg.solve(vandermonde,target)
        target[0]=0.;target[1]=1./scale
        derivative=np.linalg.solve(vandermonde,target)
        nr.extend([row]*(degree+1));nc.extend(ids%n if periodic else ids)
        nv.extend(value);dv.extend(derivative)
        # d(sum N_i f_i) = sum E_j (f_{j+1}-f_j)/h_j.
        # Edge coefficients are averages, hence the explicit interval length.
        edge_values=-np.cumsum(derivative)[:-1]*np.diff(coords)
        er.extend([row]*degree);ec.extend(ids[:-1]%count if periodic else ids[:-1]);ev.extend(edge_values)
    N=coo_matrix((nv,(nr,nc)),shape=(len(points),n)).tocsr()
    dN=coo_matrix((dv,(nr,nc)),shape=N.shape).tocsr()
    E=coo_matrix((ev,(er,ec)),shape=(len(points),count)).tocsr()
    rows=np.repeat(np.arange(count),2)
    cols=np.column_stack((np.arange(count),(np.arange(count)+1)%n)).ravel()
    data=np.column_stack((-1/np.diff(edges),1/np.diff(edges))).ravel()
    D=coo_matrix((data,(rows,cols)),shape=(count,n)).tocsr()
    return N,dN,E,D
