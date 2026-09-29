"""Continuum moving-density manufactured induction and graph-pressure targets."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import atomic_mass,electron_mass as ME,elementary_charge as E,mu_0 as MU0


@dataclass(frozen=True)
class VariableCurrentReference:
    radius: float=.035
    half_length: float=.12
    number_density: float=1e18
    contrast: float=.8
    modulation: float=.25
    tau: float=1e-5
    speed: float=200.
    axial_bias: float=.005
    amplitude: float=2e-4
    toroidal: float=.02
    ion_mass: float=39.948*atomic_mass

    @property
    def reduced_mass(self):return self.ion_mass*ME/(self.ion_mass+ME)

    def fields(self,r,z):
        r,z=np.broadcast_arrays(r,z);X=.5*self.radius**2;y=.5*r*r/X
        k=np.pi/self.half_length;sn=np.sin(k*z);cs=np.cos(k*z)
        a=self.amplitude;G=self.toroidal;eps=self.modulation
        f=(1-y)**3;h=(1-y)**2;F=1-6*y+9*y*y-4*y**3
        c=1+eps*cs;s=1+eps*sn
        density=self.number_density*(1+.3*y+self.contrast*(.6+.2*y)*sn+.1*y*cs)
        dn_r=self.number_density*(.3+.2*self.contrast*sn+.1*cs)*r/X
        dn_z=self.number_density*k*(self.contrast*(.6+.2*y)*cs-.1*y*sn)
        particle_r=self.number_density*self.speed*r/self.radius*h*cs
        particle_z=self.number_density*self.speed*(1-y)*sn
        nt=-self.number_density*self.speed*(2/self.radius*(1-4*y+3*y*y)+k*(1-y))*cs
        b=np.stack((.5*a*r*f*eps*k*sn,r*G*h*s,a*F*c),axis=-1)
        j=np.stack((-r*G*h*eps*k*cs,
                    a*r*(.5*f*eps*k*k*cs-(-6+18*y-12*y*y)*c/X),
                    2*G*(1-4*y+3*y*y)*s),axis=-1)/MU0
        electric=np.stack((np.zeros_like(y),-.5*a*r*f*c/self.tau,
                           -G*X*f*s/(3*self.tau)),axis=-1)
        coefficient=self.reduced_mass/(density*E**2)
        other=electric-coefficient[...,None]*j/self.tau+(coefficient*nt/density)[...,None]*j
        full=b.copy();full[...,2]+=self.axial_bias
        return dict(n=density,nt=nt,gradn=np.stack((dn_r,np.zeros_like(y),dn_z),axis=-1),
            particle=np.stack((particle_r,np.zeros_like(y),particle_z),axis=-1),
            psi=a*.5*r*r*f*c,g=G*h*s,magnetic=b,current=j,magnetic_rate=b/self.tau,
            current_rate=j/self.tau,electric=electric,other_forcing=other,
            current_energy=.5*coefficient*np.sum(j*j,axis=-1),
            current_energy_rate=coefficient*np.sum(j*j,axis=-1)*(1/self.tau-.5*nt/density),
            density_energy_rate=-.5*coefficient*np.sum(j*j,axis=-1)*nt/density,
            magnetic_energy_rate=np.sum(full*b/self.tau,axis=-1)/MU0)

    def thermal(self,r,z,*,gamma=5/3,density_scale=1e20):
        f=self.fields(r,z);r,z=np.broadcast_arrays(r,z);y=(r/self.radius)**2
        k=np.pi/self.half_length;sn=np.sin(k*z);cs=np.cos(k*z)
        p=10*(1+.2*y+.15*cs+.07*y*sn)
        gradp=np.stack((10*(.2+.07*sn)*2*r/self.radius**2,np.zeros_like(r),
                        10*k*(-.15*sn+.07*y*cs)),axis=-1)
        a=(f['n']/density_scale)**(gamma-1);u=p/a;w=u/f['n']
        gradw=w[...,None]*(gradp/p[...,None]-gamma*f['gradn']/f['n'][...,None])
        return dict(pressure=p,entropy=u,electric=-gradp/(E*f['n'][...,None]),
                    entropy_rate=np.sum(f['current']*gradw,axis=-1)/E)
