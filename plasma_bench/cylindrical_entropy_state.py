"""True angular-momentum state and physical-energy variables at a regular axis."""
import numpy as np
from scipy.constants import elementary_charge
from .two_entropy_flux import TwoEntropyFlux


class CylindricalEntropyState:
    """Point U=(n,n vr,n r vtheta,n vz,n se,n si)."""
    def __init__(self,total_mass,gammas):
        self.thermal=TwoEntropyFlux(total_mass,gammas);self.mass=total_mass

    def fields(self,state,radius):
        r=np.asarray(radius);u=np.asarray(state,float)
        if np.any(r<=0) or not np.isfinite(r).all():raise ValueError('Positive finite volume radii required')
        physical=u.copy();physical[...,2]/=r;f=self.thermal.fields(physical)
        gradient=f['energy_gradient'].copy();gradient[...,2]/=r
        if not np.isfinite(gradient).all():raise ValueError('Physical cylindrical energy gradient exceeds finite range')
        return dict(f,physical_state=physical,energy_gradient=gradient,omega=f['velocity'][...,1]/r)

    def from_energy_variables(self,gradient,radius):
        """Invert w=H_U at volume or boundary r, including exactly r=0."""
        w=np.asarray(gradient,float);r=np.asarray(radius);h=w[...,4:];g=self.thermal.gammas;mass=self.mass
        if not np.isfinite(w).all() or np.any(h<=0) or np.any(r<0):raise ValueError('Positive finite thermal energy-variable trace required')
        omega=w[...,2]/mass;velocity=np.stack((w[...,1]/mass,r*omega,w[...,3]/mass),axis=-1)
        log_temperature=np.log((g-1)*h/elementary_charge)
        log_density=(w[...,0]+.5*mass*np.sum(velocity*velocity,axis=-1)-np.sum(g*h,axis=-1)
                     +np.sum(h*log_temperature,axis=-1))/np.sum((g-1)*h,axis=-1)
        with np.errstate(over='ignore',under='ignore',invalid='ignore'):
            n=self.thermal.n0*np.exp(log_density);s=log_temperature-(g-1)*log_density[...,None]
            physical=np.concatenate((n[...,None],n[...,None]*velocity,n[...,None]*s),axis=-1)
        f=self.thermal.fields(physical)
        return physical,dict(f,omega=omega,energy_gradient=w)

    def pair(self,left,right,left_radius,right_radius,direction):
        """Angular flux is bar(r) Gamma bar(vtheta); symmetric in both states."""
        flux,_=self.thermal.pair(left,right,direction)
        flux[...,2]*=.5*(left_radius+right_radius)
        return flux
