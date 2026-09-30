"""Independent regular field comparison on a common cylindrical physical grid."""
import numpy as np
from scipy.signal import resample
from .cylindrical_entropy_reference import physical


def sample(driver,state,radius,axial_points):
    """Root-product r² reconstruction at arbitrary positive physical radii."""
    r=np.asarray(radius,float)
    if r.ndim!=1 or np.any(r<=0) or np.any(r>driver.radius):raise ValueError('Positive interior radial sample vector required')
    element=np.minimum((r/driver.dr).astype(int),driver.elements-1)
    x=driver.r[element,:,0]**2;y=r*r;basis=np.ones((len(r),driver.degree+1))
    for i in range(driver.degree+1):
        for j in range(driver.degree+1):
            if i!=j:basis[:,i]*=(y-x[:,j])/(x[:,i]-x[:,j])
    regular=state.copy();regular[...,1]/=driver.r;regular[...,2]/=driver.r**2
    u=np.einsum('aq,aqzk->azk',basis,regular[element]);u=resample(u,axial_points,axis=1)
    u[...,1]*=r[:,None];u[...,2]*=r[:,None]**2
    if not np.isfinite(u).all() or np.min(u[...,0])<=0:raise ValueError('Invalid regular continuous state')
    return u


def grid(radius,length,elements,points,order=24):
    q,w=np.polynomial.legendre.leggauss(order);dr=radius/elements
    r=(np.arange(elements)[:,None]*dr+.5*dr*(1+q)).ravel()
    z=np.arange(points)*length/points-.5*length
    volume=2*np.pi*.5*dr*np.tile(w,elements)*r*length/points
    return r,z,np.broadcast_to(volume[:,None],(len(r),points))


def differences(u,v,radius,volume,reference,axis_radius):
    a=physical(u,radius[:,None],reference);b=physical(v,radius[:,None],reference)
    if not all(np.isfinite(x).all() for x in (*a.values(),*b.values())) or min(np.min(a['pressure']),np.min(b['pressure']))<=0:raise ValueError('Invalid sampled physical pressure')
    result={}
    axis=radius<axis_radius
    for key in ('n','velocity','pressure','temperature','energy','entropy'):
        diff=a[key]-b[key];weight=volume if diff.ndim==2 else volume[...,None]
        den=np.sum(volume)*reference.sound**2 if key=='velocity' else np.sum(weight*b[key]**2)
        axis_den=np.sum(volume[axis])*reference.sound**2 if key=='velocity' else np.sum(weight[axis]*b[key][axis]**2)
        result[key]=dict(relative_L2=float(np.sqrt(np.sum(weight*diff**2)/den)),
                         axis_relative_L2=float(np.sqrt(np.sum(weight[axis]*diff[axis]**2)/axis_den)),
                         maximum_absolute=float(np.max(abs(diff))))
    result['continuous_energy_difference']=float((np.sum(volume*a['energy'])-np.sum(volume*b['energy']))/np.sum(volume*b['energy']))
    return result
