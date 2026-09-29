"""Independent exact-mass linear compressible two-species cylindrical waves."""
from dataclasses import dataclass
from functools import cached_property
import numpy as np
from scipy.constants import elementary_charge as E,electron_mass as ME,mu_0 as MU0
from scipy.linalg import eigh
from scipy.special import j0,j1,jvp
from .compressible_hall_reference import CompressibleHallWave,generator


@dataclass(frozen=True)
class FiniteMassWave(CompressibleHallWave):
    electron_fraction: float = 5/6
    gamma_i: float = 5/3
    electron_mass: float = ME

    def __post_init__(self):
        super().__post_init__()
        if not np.isfinite([self.electron_fraction,self.gamma_i,self.electron_mass]).all() or not 0<self.electron_fraction<1 or self.gamma_i<=1 or self.electron_mass<0 or not self.hall:
            raise ValueError('Positive species pressures, gamma_i>1, nonnegative electron mass and full physical Hall term required')

    @property
    def density(self):return self.number_density*(self.ion_mass+self.electron_mass)

    @property
    def reduced_mass(self):return self.ion_mass*self.electron_mass/(self.ion_mass+self.electron_mass)

    @property
    def delta(self):return (self.ion_mass-self.electron_mass)/(self.ion_mass+self.electron_mass)

    @property
    def electron_pressure(self):return self.electron_fraction*self.pressure

    @property
    def ion_pressure(self):return (1-self.electron_fraction)*self.pressure

    @property
    def pressure_fraction(self):return self.gamma*self.electron_pressure/(self.gamma*self.electron_pressure+self.gamma_i*self.ion_pressure)

    @property
    def sound_speed(self):return np.sqrt((self.gamma*self.electron_pressure+self.gamma_i*self.ion_pressure)/self.density)

    @property
    def inertia(self):return 1+self.reduced_mass*self.total_k**2/(MU0*self.number_density*E**2)

    @cached_property
    def eigenpair(self):
        a=generator(self.q,self.k,self.sound_speed,self.alfven_speed/np.sqrt(self.inertia),
                    self.delta*self.hall_diffusivity/self.inertia)
        values,vectors=eigh(1j*a);index=3+self.branch if self.sense==1 else 2-self.branch
        v=vectors[:,index];anchor=np.argmax(abs(v));v=self.amplitude*v*np.exp(-1j*np.angle(v[anchor]))
        v[4:]/=np.sqrt(self.inertia);v.setflags(write=False)
        return float(values[index]),v

    def dispersion(self):
        x=self.omega**2;k2=self.k**2;total=self.total_k**2;va2=self.alfven_speed**2;cs2=self.sound_speed**2
        lhs=(self.inertia*x-k2*va2)*(self.inertia*x*x-total*(self.inertia*cs2+va2)*x+k2*total*va2*cs2)
        rhs=x*self.delta**2*self.hall_diffusivity**2*k2*total*(x-total*cs2)
        return lhs,rhs

    def fields(self,r,z,time=0.,*,complex_output=False):
        out=super().fields(r,z,time,complex_output=True);r,z=np.broadcast_arrays(r,z)
        phase=np.exp(1j*(self.k*z-self.omega*time));sigma=self.eigenpair[1][3]
        amplitude=self.density*self.sound_speed*sigma*phase
        gradient=np.stack((-self.q*amplitude*j1(self.q*r),np.zeros_like(amplitude),1j*self.k*amplitude*j0(self.q*r)),axis=-1)
        fraction=self.pressure_fraction;alpha=self.electron_mass/(self.ion_mass+self.electron_mass);beta=1-alpha
        bias=np.array([0.,0.,self.bias]);current=out['current'];drift=current/(self.number_density*E)
        electric=-np.cross(out['velocity'],bias)+self.delta*np.cross(drift,bias)
        electric+=(alpha*(1-fraction)-beta*fraction)*gradient/(self.number_density*E)
        electric-=1j*self.omega*self.reduced_mass*current/(self.number_density*E**2)
        out.update(electric=electric,electron_pressure=fraction*out['pressure'],ion_pressure=(1-fraction)*out['pressure'],
            pressure_gradient=gradient,electron_velocity=out['velocity']-beta*drift,ion_velocity=out['velocity']+alpha*drift)
        return out if complex_output else {key:value.real for key,value in out.items()}

    def derivatives(self,r,z,time=0.):
        """Analytic Bessel/Fourier differentiation, independent of amplitude equations."""
        r,z=np.broadcast_arrays(r,z);anchor=self.radius/3
        fields=self.fields(r,z,time,complex_output=True)
        sample=self.fields(np.full_like(r,anchor),z,time,complex_output=True)
        dr={}
        for key,value in fields.items():
            if value.ndim==r.ndim+1:
                dr[key]=sample[key]*np.stack((self.q*jvp(1,self.q*r)/j1(self.q*anchor),
                    self.q*jvp(1,self.q*r)/j1(self.q*anchor),-self.q*j1(self.q*r)/j0(self.q*anchor)),axis=-1)
            elif key=='psi':dr[key]=sample[key]*self.q*r*j0(self.q*r)/(anchor*j1(self.q*anchor))
            else:dr[key]=sample[key]*(-self.q*j1(self.q*r)/j0(self.q*anchor))
        return dr,{key:1j*self.k*value for key,value in fields.items()},{key:-1j*self.omega*value for key,value in fields.items()}
