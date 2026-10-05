"""Original two-species energy, entropy, particle and current fluxes."""
import numpy as np
from scipy.constants import elementary_charge as E,electron_mass as ME,mu_0 as MU0


def primitives(sample,state):
    n,v,_=sample.fluid(state);d=sample.d
    entropy=np.stack([sample.scalar(state[...,k])/n for k in (4,5)],axis=-1)
    h=E*np.exp(entropy)*(n[...,None]/d.thermal.n0)**(d.thermal.gammas-1)/(d.thermal.gammas-1)
    return n,v,entropy,h


def fluxes(sample,state,magnetic,electric,current,bias,ion_mass):
    n,v,entropy,h=primitives(sample,state)
    j=sample.electric(current);e=sample.electric(electric);b=sample.magnetic(magnetic).copy();b[...,2]+=bias
    mass=ion_mass+ME;alpha=ME/mass;beta=ion_mass/mass
    velocities=[v-beta*j/(E*n[...,None]),v+alpha*j/(E*n[...,None])]
    energy=np.cross(e,b)/MU0;entropy_flux=[]
    for k,(velocity,species_mass) in enumerate(zip(velocities,(ME,ion_mass))):
        energy+=velocity*(.5*species_mass*n*np.sum(velocity*velocity,axis=-1)+sample.d.thermal.gammas[k]*n*h[...,k])[...,None]
        entropy_flux.append(velocity*(n*entropy[...,k])[...,None])
    return np.stack((n[...,None]*v,j,*entropy_flux,energy),axis=-2)


def rates(sample,state,magnetic,current,rate,bias,ion_mass):
    n,v,_,nt,vt=sample.fluid(state,rate.fluid)
    entropy=np.stack([sample.scalar(state[...,k])/n for k in (4,5)],axis=-1)
    ct=np.stack([sample.scalar(rate.fluid[...,k]) for k in (4,5)],axis=-1)
    gammas=sample.d.thermal.gammas
    h=E*np.exp(entropy)*(n[...,None]/sample.d.thermal.n0)**(gammas-1)/(gammas-1)
    j,dj=sample.electric(current,True);jt=sample.electric(rate.current)
    b=sample.magnetic(magnetic).copy();b[...,2]+=bias;bt=sample.magnetic(rate.magnetic)
    mass=ion_mass+ME;reduced=ion_mass*ME/mass
    thermal=h*(ct+(gammas-entropy)*nt[...,None])
    power=dict(bulk_kinetic=mass*(.5*nt*np.sum(v*v,axis=-1)+n*np.sum(v*vt,axis=-1)),
        electron_internal=thermal[...,0],ion_internal=thermal[...,1],
        current=reduced/E**2*(np.sum(j*jt,axis=-1)/n-.5*np.sum(j*j,axis=-1)*nt/n**2),
        magnetic=np.sum(b*bt,axis=-1)/MU0)
    conserved=np.stack((nt,np.zeros_like(nt),ct[...,0],ct[...,1],sum(power.values())),axis=-1)
    return conserved,power,np.trace(dj,axis1=-2,axis2=-1),np.sum(dj*dj,axis=(-1,-2))
