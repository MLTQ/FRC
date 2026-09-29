"""Exact periodic thermal transport by a radially sheared compatible current."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E,mu_0 as MU0
from scipy.special import i0


@dataclass(frozen=True)
class LogEntropyReference:
    radius: float=.035
    half_length: float=.12
    number_density: float=1e18
    density_slope: float=.3
    toroidal: float=.02
    background_entropy: float=2.
    amplitude: float=.3
    gamma: float=5/3
    density_scale: float=1e20
    ion_mass: float=39.948*atomic_mass

    @property
    def axis_speed(self):return 2*self.toroidal/(MU0*E*self.number_density)

    @property
    def transit_time(self):return 2*self.half_length/abs(self.axis_speed)

    def fields(self,r,z,time=0.):
        r,z=np.broadcast_arrays(r,z);y=(r/self.radius)**2;k=np.pi/self.half_length
        n=self.number_density*(1+self.density_slope*y)
        jz=2*self.toroidal*(1-4*y+3*y*y)/MU0
        speed=jz/(E*n);phase=k*(z+speed*time)
        s=self.background_entropy+self.amplitude*np.cos(phase)
        st=-self.amplitude*k*speed*np.sin(phase)
        w=E*np.exp(s);pressure=n*(n/self.density_scale)**(self.gamma-1)*w
        return dict(n=n,g=self.toroidal*(1-y)**2,entropy=s,entropy_rate=st,
            pressure=pressure,internal=pressure/(self.gamma-1),proxy=n*w,
            current=np.stack((np.zeros_like(r),np.zeros_like(r),jz),axis=-1))

    def integrals(self):
        volume=2*np.pi*self.radius**2*self.half_length;d=self.density_slope
        particles=volume*self.number_density*(1+.5*d)
        mean_power=((1+d)**(self.gamma+1)-1)/(d*(self.gamma+1)) if d else 1.
        factor=E*np.exp(self.background_entropy)*i0(self.amplitude)
        energy=volume*self.number_density**self.gamma/self.density_scale**(self.gamma-1)*mean_power*factor/(self.gamma-1)
        return dict(particles=particles,entropy=particles*self.background_entropy,
                    proxy=particles*factor,internal=energy)
