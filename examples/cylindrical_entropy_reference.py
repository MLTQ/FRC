"""Independent cylindrical physical fields and regular continuous reconstruction."""
import numpy as np
from scipy.constants import elementary_charge as E
from scipy.signal import resample


def initial(driver,reference,time=0.):
    f=reference.fields(driver.r,driver.z,time);n=f['n']
    u=np.concatenate((n[...,None],n[...,None]*f['velocity'],
                      n[...,None]*np.stack((f['electron'],f['ion']),axis=-1)),axis=-1)
    u[...,2]*=driver.r
    return u


def physical(u,r,reference):
    """Original equation of state; supports complex-step differentiation."""
    n=u[...,0];v=u[...,1:4]/n[...,None];v=v.copy();v[...,1]/=r
    s=u[...,4:]/n[...,None];g=np.array([reference.gamma_e,reference.gamma_i])
    temperature=np.exp(s)*(n[...,None]/1e20)**(g-1);p=n[...,None]*E*temperature
    return dict(n=n,velocity=v,entropy=s,temperature=temperature,pressure=p,
                energy=.5*reference.mass*n*np.sum(v*v,axis=-1)+np.sum(p/(g-1),axis=-1))


def rates(u,du,r,reference):
    f=physical(u,r,reference);n=f['n'];v=f['velocity'];s=f['entropy']
    dm=du[...,1:4].copy();dm[...,1]/=r
    vt=(dm-v*du[...,0,None])/n[...,None]
    st=(du[...,4:]-s*du[...,0,None])/n[...,None]
    g=np.array([reference.gamma_e,reference.gamma_i])
    pt=f['pressure']*(st+g*du[...,0,None]/n[...,None])
    et=reference.mass*(.5*du[...,0]*np.sum(v*v,axis=-1)+n*np.sum(v*vt,axis=-1))+np.sum(pt/(g-1),axis=-1)
    return dict(n=du[...,0],velocity=vt,entropy=st,pressure=pt,energy=et)


def regular_sample(driver,u,radial_order=32,axial_factor=4):
    """Independent root-product interpolation in r², with regular momenta."""
    q,w=np.polynomial.legendre.leggauss(radial_order)
    r=(np.arange(driver.elements)[:,None]*driver.dr+.5*driver.dr*(1+q))[...,None]
    # Normalized local r² avoids a large monomial Vandermonde system. Each
    # product is an independent Lagrange basis, not the production SBP matrix.
    a=(np.arange(driver.elements)*driver.dr)**2
    width=((np.arange(driver.elements)+1)*driver.dr)**2-a
    x=(driver.r[...,0]**2-a[:,None])/width[:,None]
    y=(r[...,0]**2-a[:,None])/width[:,None]
    basis=np.ones((driver.elements,radial_order,driver.degree+1))
    for i in range(driver.degree+1):
        for j in range(driver.degree+1):
            if i!=j:basis[:,:,i]*=(y-x[:,j,None])/(x[:,i,None]-x[:,j,None])
    regular=u.copy();regular[...,1]/=driver.r;regular[...,2]/=driver.r**2
    sample=np.einsum('eaq,eqzk->eazk',basis,regular)
    sample=resample(sample,driver.nz*axial_factor,axis=2)
    sample[...,1]*=r;sample[...,2]*=r*r
    z=(np.arange(driver.nz*axial_factor)*driver.dz/axial_factor-.5*driver.length)[None,None,:]
    volume=np.broadcast_to(2*np.pi*(driver.dz/axial_factor)*.5*driver.dr*w[None,:,None]*r,sample.shape[:-1])
    if not np.isfinite(sample).all() or np.any(sample[...,0]<=0):raise ValueError('Regular physical reconstruction has invalid density')
    return sample,r,z,volume


def measure(driver,u,reference,time=0.,radial_order=32,axial_factor=4):
    sample,r,z,volume=regular_sample(driver,u,radial_order,axial_factor)
    f=physical(sample,r,reference);exact=reference.fields(r,z,time)
    targets=dict(n=exact['n'],velocity=exact['velocity'],pressure=np.stack((exact['pe'],exact['pi']),axis=-1),
                 temperature=np.stack((exact['Te'],exact['Ti']),axis=-1),energy=exact['energy'])
    if not all(np.isfinite(v).all() for v in f.values()) or np.any(f['pressure']<=0):raise ValueError('Regular reconstructed thermodynamic fields are invalid')
    metrics={}
    for key,wanted in targets.items():
        weight=volume if wanted.ndim==3 else volume[...,None]
        denominator=np.sum(volume)*reference.sound**2 if key=='velocity' else np.sum(weight*wanted**2)
        name='velocity_over_sound' if key=='velocity' else key
        metrics[name]=float(np.sqrt(np.sum(weight*(f[key]-wanted)**2)/denominator))
    continuous=float(np.sum(volume*f['energy']));nodal=float(np.sum(driver.volume*physical(u,driver.r,reference)['energy']))
    metrics.update(continuous_energy=continuous,nodal_energy=nodal,quadrature_relative_difference=(nodal-continuous)/continuous,
                   minimum_density=float(np.min(f['n'])),minimum_pressure=float(np.min(f['pressure'])))
    return metrics
