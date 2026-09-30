"""Independent original axial Euler equations and continuous physical sampling."""
import numpy as np
from scipy.constants import elementary_charge as E
from scipy.signal import resample
from .log_fluid_reference import LogFluidReference


def state(reference,z,time=0.,manufactured=False):
    if not manufactured:
        f=reference.fields(np.zeros_like(z),z,time);n=f['n'];v=f['velocity'];s=np.stack((f['electron'],f['ion']),axis=-1)
    else:
        phase=np.pi*z/reference.half_length;n=reference.number_density*(1+.3*np.cos(phase));c=reference.sound
        v=c*np.stack((.1*np.sin(2*phase),.15*np.cos(2*phase),.2+.1*np.sin(phase)),axis=-1)
        p=reference.pressure*np.stack((.6*(1+.2*np.cos(phase)),.4*(1+.15*np.sin(2*phase))),axis=-1)
        gamma=np.array([reference.gamma_e,reference.gamma_i]);s=np.log(p/(n[...,None]*E*(n[...,None]/1e20)**(gamma-1)))
    return np.concatenate((n[...,None],n[...,None]*v,n[...,None]*s),axis=-1)


def physical(u,reference):
    n=u[...,0];v=u[...,1:4]/n[...,None];s=u[...,4:]/n[...,None];gamma=np.array([reference.gamma_e,reference.gamma_i])
    temperature=np.exp(s)*(n[...,None]/1e20)**(gamma-1);p=n[...,None]*E*temperature
    return dict(n=n,velocity=v,entropy=s,temperature=temperature,pressure=p,
                energy=.5*reference.mass*n*np.sum(v*v,axis=-1)+np.sum(p/(gamma-1),axis=-1))


def parent_rate(reference,z,manufactured=False):
    u=state(reference,z,manufactured=manufactured);f=physical(u,reference);h=1e-28
    derivative=np.imag(state(reference,z+1j*h,manufactured=manufactured))/h
    differentiated=physical(state(reference,z+1j*h,manufactured=manufactured),reference)
    dz={k:np.imag(v)/h for k,v in differentiated.items()}
    result=-f['velocity'][...,2,None]*derivative-u*dz['velocity'][...,2,None]
    result[...,3]-=np.sum(dz['pressure'],axis=-1)/reference.mass
    energy_rate=-(dz['velocity'][...,2]*(f['energy']+np.sum(f['pressure'],axis=-1))+
                  f['velocity'][...,2]*(dz['energy']+np.sum(dz['pressure'],axis=-1)))
    return result,energy_rate


def energy_rate(u,rate,reference):
    f=physical(u,reference);n=f['n'];v=f['velocity'];gamma=np.array([reference.gamma_e,reference.gamma_i])
    nt=rate[...,0];st=(rate[...,4:]-f['entropy']*nt[...,None])/n[...,None]
    kinetic=reference.mass*(np.sum(v*rate[...,1:4],axis=-1)-.5*np.sum(v*v,axis=-1)*nt)
    return kinetic+np.sum(f['pressure']*(st+gamma*nt[...,None]/n[...,None])/(gamma-1),axis=-1)


def measure(u,reference,time,oversample=8):
    size=len(u)*oversample;sample=resample(u,size,axis=0);z=np.arange(size)*(2*reference.half_length)/size-reference.half_length
    if np.min(sample[:,0])<=0:raise ValueError('Independent continuous density is nonpositive')
    f=physical(sample,reference);exact=physical(state(reference,z,time),reference)
    metrics={key:float(np.linalg.norm(f[key]-exact[key])/np.linalg.norm(exact[key]))
             for key in ('n','velocity','pressure','temperature','energy')}
    metrics['density_signal']=float(np.linalg.norm(f['n']-exact['n'])/max(np.linalg.norm(exact['n']-reference.number_density),1e-300))
    metrics['continuous_energy']=float(np.mean(f['energy'])*2*reference.half_length)
    metrics['nodal_energy']=float(np.mean(physical(u,reference)['energy'])*2*reference.half_length)
    metrics['quadrature_relative_difference']=(metrics['nodal_energy']-metrics['continuous_energy'])/metrics['continuous_energy']
    return metrics
