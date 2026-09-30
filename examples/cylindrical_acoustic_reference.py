"""Small-amplitude Bessel standing wave of the original cylindrical Euler PDE."""
from dataclasses import dataclass
import numpy as np
from scipy.special import j0,j1,jn_zeros
from scipy.constants import elementary_charge as E
from .log_fluid_reference import LogFluidReference


@dataclass(frozen=True)
class CylindricalAcousticReference(LogFluidReference):
    kind: str='radial_acoustic'
    gamma_i: float=5/3
    amplitude: float=1e-4

    @property
    def duration(self):return 2*np.pi*self.radius/(jn_zeros(1,1)[0]*self.sound)

    def fields(self,r,z,time=0.):
        r,z=np.broadcast_arrays(r,z);k=jn_zeros(1,1)[0]/self.radius
        n=self.number_density*(1+self.amplitude*j0(k*r)*np.cos(k*self.sound*time))
        vr=self.amplitude*self.sound*j1(k*r)*np.sin(k*self.sound*time)
        velocity=np.stack((vr,np.zeros_like(r),np.zeros_like(r)),axis=-1)
        pe=.6*self.pressure*(n/self.number_density)**self.gamma_e;pi=.4*self.pressure*(n/self.number_density)**self.gamma_i
        s=np.log(np.stack((pe,pi),axis=-1)/(n[...,None]*E*(n[...,None]/1e20)**(np.array([self.gamma_e,self.gamma_i])-1)))
        return dict(n=n,velocity=velocity,electron=s[...,0],ion=s[...,1],pe=pe,pi=pi,Te=pe/(n*E),Ti=pi/(n*E),
                    energy=.5*self.mass*n*vr*vr+pe/(self.gamma_e-1)+pi/(self.gamma_i-1))
