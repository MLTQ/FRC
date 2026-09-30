"""Joint smooth-fluid flux conserving both entropy densities and physical energy."""
import numpy as np
from scipy.constants import elementary_charge
from .thermal_means import logarithmic_mean


class TwoEntropyFlux:
    """Point state (n, n*vx, n*vy, n*vz, n*se, n*si); not cell means."""
    def __init__(self,total_mass,gammas=(5/3,1.4),reference_density=1e20):
        self.mass=float(total_mass);self.gammas=np.asarray(gammas,float);self.n0=float(reference_density)
        if self.gammas.shape!=(2,) or not np.isfinite(self.gammas).all() or np.any(self.gammas<=1):
            raise ValueError('Two finite species gammas greater than one required')
        if not np.isfinite([self.mass,self.n0]).all() or min(self.mass,self.n0)<=0:
            raise ValueError('Positive finite physical mass and reference density required')

    def fields(self,state):
        u=np.asarray(state,float)
        if u.shape[-1:]!=(6,) or not np.isfinite(u).all() or np.any(u[...,0]<=0):
            raise ValueError('Finite six-component state with positive density required')
        n=u[...,0];v=u[...,1:4]/n[...,None];s=u[...,4:6]/n[...,None]
        with np.errstate(over='ignore',under='ignore',invalid='ignore'):
            temperature=np.exp(s+(self.gammas-1)*np.log(n[...,None]/self.n0))
            h=elementary_charge*temperature/(self.gammas-1)
            pressure=(self.gammas-1)*n[...,None]*h
            kinetic=.5*self.mass*n*np.sum(v*v,axis=-1);energy=kinetic+n*np.sum(h,axis=-1)
        if not np.isfinite(energy).all() or not np.isfinite(v).all() or not np.isfinite(h).all() or np.any(h<=0) or np.any(pressure<=0):
            raise ValueError('Physical entropy/kinetic state exceeds positive finite range')
        gradient=np.concatenate(((-.5*self.mass*np.sum(v*v,axis=-1)+np.sum((self.gammas-s)*h,axis=-1))[...,None],
                                 self.mass*v,h),axis=-1)
        return dict(n=n,velocity=v,entropy=s,specific_internal=h,temperature=temperature,pressure=pressure,
                    kinetic=kinetic,energy=energy,energy_gradient=gradient)

    def pair(self,left,right,direction=2):
        if direction not in (0,1,2):raise ValueError('Cartesian flux direction 0,1,2 required')
        l=self.fields(left);r=self.fields(right);nl=logarithmic_mean(l['n'],r['n']);na=.5*l['n']+.5*r['n']
        va=.5*l['velocity']+.5*r['velocity'];speed=va[...,direction];mass=nl*speed
        ha=.5*l['specific_internal']+.5*r['specific_internal'];hl=logarithmic_mean(l['specific_internal'],r['specific_internal'])
        sa=.5*l['entropy']+.5*r['entropy'];pa=.5*np.sum(l['pressure'],axis=-1)+.5*np.sum(r['pressure'],axis=-1)
        momentum=mass[...,None]*va;momentum[...,direction]+=pa/self.mass
        entropy=speed[...,None]*(nl[...,None]*sa+(self.gammas-1)*(na-nl)[...,None]+nl[...,None]*(ha/hl-1))
        flux=np.concatenate((mass[...,None],momentum,entropy),axis=-1)
        potential_l=np.sum(l['pressure'],axis=-1)*l['velocity'][...,direction]
        potential_r=np.sum(r['pressure'],axis=-1)*r['velocity'][...,direction]
        energy_flux=np.sum((.5*l['energy_gradient']+.5*r['energy_gradient'])*flux,axis=-1)-.5*(potential_l+potential_r)
        if not np.isfinite(flux).all() or not np.isfinite(energy_flux).all():raise ValueError('Physical flux exceeds finite range')
        return flux,energy_flux
