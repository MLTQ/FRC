"""Exact collisionless sphere response to an initially zero uniform bias."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0, electron_mass as ME, elementary_charge as E
from scipy.special import ive


@dataclass(frozen=True)
class InertialScreening:
    radius: float = 1.
    electron_length: float = .2
    bias: float = 1e-8

    def __post_init__(self):
        if not np.isfinite(self.radius+self.electron_length+self.bias) or min(self.radius,self.electron_length)<=0:
            raise ValueError("Finite positive sphere/electron lengths and finite bias required")

    @property
    def electron_density(self):
        return ME/(MU0*E**2*self.electron_length**2)

    def radial(self,s):
        s=np.asarray(s,float)
        if not np.all(np.isfinite(s)) or np.any(s<0):raise ValueError("Finite nonnegative radius required")
        xb=self.radius/self.electron_length;x=np.minimum(s,self.radius)/self.electron_length
        small=x<1e-3;safe=np.where(small,1.,x)
        inv_i0=2*xb*np.exp(-xb)/(-np.expm1(-2*xb))
        scale=np.sqrt(xb/safe)*np.exp(safe-xb)/ive(.5,xb)
        f=scale*ive(1.5,safe)/safe
        derivative=(scale*ive(.5,safe)-3*f)/safe
        f=np.where(small,inv_i0*(1/3+x*x/30+x**4/840+x**6/45360),f)
        derivative=np.where(small,inv_i0*(x/15+x**3/210+x**5/7560),derivative)
        interior=1.5*self.bias*f;inside_derivative=1.5*self.bias/self.electron_length*derivative
        boundary=1.5*self.bias*ive(1.5,xb)/(xb*ive(.5,xb))
        exterior_self=(boundary-.5*self.bias)*(self.radius/np.maximum(s,self.radius))**3
        return (np.where(s<=self.radius,interior,.5*self.bias+exterior_self),
                np.where(s<=self.radius,inside_derivative,-3*exterior_self/np.maximum(s,self.radius)))

    def fields(self,r,z):
        r,z=np.broadcast_arrays(np.asarray(r,float),np.asarray(z,float));s=np.hypot(r,z)
        a,derivative=self.radial(s)
        slope=np.divide(derivative,s,out=np.zeros_like(s),where=s!=0)
        return {"a":a,"Br":-r*z*slope,"Bz":2*a+r*r*slope,
                "Jtheta":np.where(s<=self.radius,-r*a/(MU0*self.electron_length**2),0.)}
