"""Independent conducting-sphere eigenmode with exact vacuum, current and heat."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.special import spherical_jn


@dataclass(frozen=True)
class SphericalDecay:
    radius: float = 1.
    diffusivity: float = .1
    central_a: float = float(np.sqrt(MU0))
    mode: int = 1

    def __post_init__(self):
        if not np.isfinite(self.radius+self.diffusivity+self.central_a) or self.radius <= 0 or self.diffusivity < 0:
            raise ValueError("Finite positive radius and nonnegative diffusivity required")
        if isinstance(self.mode,bool) or not isinstance(self.mode,int) or self.mode < 1:
            raise ValueError("Positive integer radial mode required")

    @property
    def wavenumber(self):
        return self.mode*np.pi/self.radius

    def radial(self,s,time=0.):
        if not np.isfinite(time) or time < 0:
            raise ValueError("Finite nonnegative reference time required")
        s=np.asarray(s,float)
        if np.any(s < 0) or not np.all(np.isfinite(s)):
            raise ValueError("Finite nonnegative spherical radius required")
        k=self.wavenumber;x=k*s
        regular=np.abs(x)<1e-3;safe=np.where(regular,1.,x)
        f=np.where(regular,1-x*x/10+x**4/280-x**6/15120,3*spherical_jn(1,safe)/safe)
        df=np.where(regular,-x/5+x**3/70-x**5/2520,
                    3*(spherical_jn(1,safe,True)/safe-spherical_jn(1,safe)/safe**2))
        decay=np.exp(-self.diffusivity*k*k*time)
        boundary=3*self.central_a*(-1)**(self.mode+1)/(self.mode*np.pi)**2
        outside=boundary*(self.radius/np.maximum(s,self.radius))**3
        a=np.where(s <= self.radius,self.central_a*f,outside)*decay
        derivative=np.where(s <= self.radius,self.central_a*k*df,-3*outside/np.maximum(s,self.radius))*decay
        return a,derivative

    def fields(self,r,z,time=0.):
        r,z=np.broadcast_arrays(np.asarray(r,float),np.asarray(z,float))
        s=np.hypot(r,z);a,derivative=self.radial(s,time)
        slope=np.divide(derivative,s,out=np.full_like(s,-self.central_a*self.wavenumber**2/5*
                         np.exp(-self.diffusivity*self.wavenumber**2*time)),where=s!=0)
        current=np.where(s <= self.radius,self.wavenumber**2/MU0*r*a,0.)
        return {"a":a,"psi":r*r*a,"Br":-r*z*slope,"Bz":2*a+r*r*slope,
                "Jtheta":current,"Etheta":self.diffusivity*self.wavenumber**2*r*a}

    def heat_density(self,r,z,time):
        if not np.isfinite(time) or time < 0:
            raise ValueError("Finite nonnegative reference time required")
        current=self.fields(r,z,0.)["Jtheta"]
        fraction=-np.expm1(-2*self.diffusivity*self.wavenumber**2*time)
        return MU0*current**2/(2*self.wavenumber**2)*fraction

    def integrals(self,time=0.):
        if not np.isfinite(time) or time < 0:
            raise ValueError("Finite nonnegative reference time required")
        initial=6*np.pi*self.central_a**2*self.radius/(MU0*self.wavenumber**2)
        rate=self.diffusivity*self.wavenumber**2;decay=np.exp(-2*rate*time)
        return {"induced_energy_j":float(initial*decay),"ohmic_power_w":float(2*rate*initial*decay),
                "ohmic_heat_j":float(-initial*np.expm1(-2*rate*time))}
