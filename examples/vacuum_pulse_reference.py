"""Independent cylindrical diffusion limits for a plasma surrounded by vacuum."""
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.special import j0, j1, jn_zeros, kv
from scipy.optimize import brentq
from .flux_diffusion_reference import FluxDiffusionReference


class VacuumPulseReference(FluxDiffusionReference):
    def __init__(self, radial_points, **options):
        super().__init__(radial_points,**options)
        roots=jn_zeros(0,len(self.q)); self.q=roots/self.radius
        self.constant_coefficients=2/(roots*j1(roots))
        self.rates=self.diffusivity*self.q**2/(1+self.de2*self.q**2)
        self.basis_b=j0(self.q[:,None]*self.r.ravel()[None,:])
        self.basis_j=self.q[:,None]*j1(self.q[:,None]*self.r.ravel()[None,:])/MU0

    def evaluate(self, when):
        t=float(when)
        source=self.amplitude*np.sin(np.pi*t/self.period)**2 if 0<t<self.period else 0.
        a=float(self.convolution(self.diffusivity/self.de2,t)/(4*self.de2))
        modal=-self.constant_coefficients*self.convolution(self.rates,t)/(1+self.de2*self.q**2)
        residual=modal+4*a*self.constant_coefficients/self.q**2
        bz=source+a*(self.r**2-self.radius**2)+(residual@self.basis_b).reshape(self.r.shape)
        jtheta=-2*a*self.r/MU0+(residual@self.basis_j).reshape(self.r.shape)
        return bz,jtheta


class VacuumDecayMode:
    def __init__(self, reference, *, length, amplitude):
        self.radius=reference.radius; self.k=2*np.pi/length; self.amplitude=amplitude
        x=self.k*self.radius; robin=x*kv(0,x)/kv(1,x)
        root=brentq(lambda y:y*j0(y)+robin*j1(y),jn_zeros(0,1)[0],jn_zeros(1,1)[0])
        self.q=root/self.radius
        eigenvalue=self.q**2+self.k**2
        self.rate=reference.diffusivity*eigenvalue/(1+reference.de2*eigenvalue)

    def flux(self, r, z):
        return self.amplitude/self.q*np.asarray(r)*j1(self.q*np.asarray(r))*np.cos(self.k*np.asarray(z))

    def plasma(self, r, z, when):
        a=self.amplitude*np.exp(-self.rate*when); r=np.asarray(r); z=np.asarray(z)
        br=a*self.k/self.q*j1(self.q*r)*np.sin(self.k*z)
        bz=a*j0(self.q*r)*np.cos(self.k*z)
        jtheta=a*(self.q**2+self.k**2)/(MU0*self.q)*j1(self.q*r)*np.cos(self.k*z)
        return br,bz,jtheta

    def vacuum(self, r, z, when):
        r=np.asarray(r); z=np.asarray(z); k=self.k; R=self.radius
        trace=self.amplitude/self.q*R*j1(self.q*R)*np.exp(-self.rate*when)
        br=k*trace*kv(1,k*r)/(R*kv(1,k*R))*np.sin(k*z)
        bz=-k*trace*kv(0,k*r)/(R*kv(1,k*R))*np.cos(k*z)
        return br,bz
