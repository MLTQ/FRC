"""Regular physical fields from the active nonsingular-coordinate DOFs."""
import numpy as np
from scipy.signal import resample


def basis(nodes,points):
    matrix=np.ones((len(points),len(nodes)))
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if i!=j:matrix[:,i]*=(points-nodes[j])/(nodes[i]-nodes[j])
    return matrix


def sample(driver,state,radius,axial_points):
    r=np.asarray(radius,float);element=np.minimum((r*driver.elements/driver.radius).astype(int),driver.elements-1)
    if r.ndim!=1 or np.min(r)<0 or np.max(r)>driver.radius:raise ValueError('Physical radial samples outside cylinder')
    result=np.empty((len(r),driver.nz,6));velocity=np.empty((len(r),driver.nz,2))
    for e in range(driver.elements):
        select=np.flatnonzero(element==e)
        if not len(select):continue
        x=driver.x[e,:,0];target=.5*r[select]**2;interpolation=basis(x,target)
        result[select]=np.einsum('aq,qzk->azk',interpolation,state[e])
        begin=1 if e==0 else 0
        v=state[e,begin:,:,1:3]/state[e,begin:,:,0,None]
        velocity[select]=np.einsum('aq,qzk->azk',basis(x[begin:],target),v)
    # Interpolate density/mz/C linearly, regular transverse velocities in their
    # declared polynomial spaces, then form their actual physical momenta.
    result=resample(result,axial_points,axis=1);velocity=resample(velocity,axial_points,axis=1)
    n=result[...,0];result[...,1]=n*r[:,None]*velocity[...,0];result[...,2]=n*r[:,None]**2*velocity[...,1]
    if not np.isfinite(result).all() or np.any(n<=0):raise ValueError('Invalid regular physical reconstruction')
    return result
