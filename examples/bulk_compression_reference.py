"""Independent exact periodic adiabatic compression and translation targets."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import atomic_mass,elementary_charge as E


@dataclass(frozen=True)
class BulkCompressionReference:
    radius: float=.035
    half_length: float=.12
    number_density: float=1e18
    ion_mass: float=39.948*atomic_mass
    speed: float=2000.
    density_contrast: float=.25
    gamma_e: float=5/3
    gamma_i: float=1.4
    kind: str='compression'

    @property
    def duration(self):
        return np.log(2)*self.half_length/(np.pi*abs(self.speed)) if self.kind=='compression' else 2*self.half_length/abs(self.speed)

    def velocity(self,z):
        return self.speed*np.sin(np.pi*z/self.half_length) if self.kind=='compression' else np.zeros_like(z)+self.speed

    def fields(self,r,z,time=0.):
        r,z=np.broadcast_arrays(r,z);k=np.pi/self.half_length;y=(r/self.radius)**2
        if self.kind=='compression':
            factor=np.exp(-k*self.speed*time);c=np.cos(.5*k*z);s=np.sin(.5*k*z)
            origin=2/k*np.arctan(factor*np.tan(.5*k*z));jac=factor/(c*c+factor**2*s*s)
        elif self.kind=='translation':origin=z-self.speed*time;jac=np.ones_like(z)
        else:raise ValueError('Declared compression or translation reference required')
        n=self.number_density*(1+.3*y)*(1+self.density_contrast*np.cos(k*origin))*jac
        electron=3+.1*y+.3*np.sin(k*origin)
        ion=2.5-.15*y+.2*np.cos(k*origin+.23)
        result=dict(n=n,electron=electron,ion=ion,velocity=self.velocity(z))
        for label,entropy,gamma in (('electron',electron,self.gamma_e),('ion',ion,self.gamma_i)):
            temperature=np.exp(entropy)*(n/1e20)**(gamma-1)
            result[label+'_temperature']=temperature
            result[label+'_pressure']=n*E*temperature
            result[label+'_internal']=n*E*temperature/(gamma-1)
        return result

    def momentum(self,model,time):
        # Uses the evolved numerical density only. Analytic n(r,z,t) is never
        # passed into an evolution stage.
        if self.kind=='translation':axial=model.density_map.cell_means*self.speed
        else:axial=model.average(model.density*self.velocity(model.geometry.z))
        return np.stack((np.zeros_like(axial),axial),axis=-1)
