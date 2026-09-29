"""Exact singly charged quasineutral species-to-bulk energy and momentum moments."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import electron_mass as ME, elementary_charge as E, mu_0 as MU0


@dataclass(frozen=True)
class SpeciesMasses:
    ion: float
    electron: float = ME

    def __post_init__(self):
        if not np.isfinite([self.ion,self.electron]).all() or self.ion<=0 or self.electron<0:
            raise ValueError('Positive ion mass and nonnegative physical electron mass required')

    @property
    def total(self):return self.ion+self.electron

    @property
    def alpha(self):return self.electron/self.total

    @property
    def beta(self):return self.ion/self.total

    @property
    def delta(self):return self.beta-self.alpha

    @property
    def reduced(self):return self.ion*self.alpha


def moments(masses,n,velocity,current,magnetic,pe,pi,electric,*,gamma_e=5/3,gamma_i=5/3):
    """Local SI moments and flux; vectors use the physical orthonormal basis."""
    n,pe,pi=np.broadcast_arrays(n,pe,pi)
    if not np.isfinite([gamma_e,gamma_i]).all() or min(gamma_e,gamma_i)<=1:
        raise ValueError('Species gammas must exceed one')
    if any(not np.isfinite(x).all() for x in (n,pe,pi,velocity,current,magnetic,electric)) or np.any(n<=0) or np.any(pe<=0) or np.any(pi<=0):
        raise ValueError('Finite fields and strictly positive density/species pressures required')
    v,j,b,ef=(np.broadcast_to(x,n.shape+(3,)) for x in (velocity,current,magnetic,electric))
    w=j/(n[...,None]*E);rho=masses.total*n
    bulk=.5*rho*np.sum(v*v,axis=-1);drift=.5*masses.reduced*n*np.sum(w*w,axis=-1)
    he=gamma_e*pe/(gamma_e-1);hi=gamma_i*pi/(gamma_i-1)
    stress=masses.reduced*n[...,None,None]*w[..., :,None]*w[...,None,:]
    flux=(bulk+drift+he+hi)[...,None]*v+np.cross(ef,b)/MU0
    flux+=w*(masses.reduced*n*np.sum(v*w,axis=-1)-masses.delta*drift+masses.alpha*hi-masses.beta*he)[...,None]
    return dict(drift=w,ion_velocity=v+masses.alpha*w,electron_velocity=v-masses.beta*w,
        bulk_kinetic=bulk,current_kinetic=drift,drift_stress=stress,
        total_energy=bulk+drift+pe/(gamma_e-1)+pi/(gamma_i-1)+np.sum(b*b,axis=-1)/(2*MU0),
        total_energy_flux=flux)


def reduced_ohm(masses,n,v,w,b,grad_pe,grad_pi,w_rate,adv_v_w,adv_w_v,adv_w_w):
    """Full continuum Ohm field, given actual drift time/convective derivatives."""
    return (-np.cross(v,b)+masses.delta*np.cross(w,b)
        +(masses.alpha*grad_pi-masses.beta*grad_pe)/(n[...,None]*E)
        +masses.reduced/E*(w_rate+adv_v_w+adv_w_v-masses.delta*adv_w_w))
