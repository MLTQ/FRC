"""Independent analytic axisymmetric moving two-species manufactured fields."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import atomic_mass, electron_mass as ME, elementary_charge as E, mu_0 as MU0


@dataclass(frozen=True)
class MovingSpeciesReference:
    radius: float = .035
    length: float = .24
    density: float = 1e18
    ion_mass: float = 39.948*atomic_mass
    electron_mass: float = ME
    frequency: float = 1e7
    bias: float = .005
    poloidal: float = .0006
    toroidal: float = .045
    electron_ev: float = 5.
    ion_ev: float = 1.
    gamma_e: float = 5/3
    gamma_i: float = 1.4

    def fields(self,r,z,t):
        r,z,t=np.broadcast_arrays(r,z,t);s=(r/self.radius)**2;k=2*np.pi/self.length
        phase=k*z-.7*self.frequency*t;other=k*z+.2-.9*self.frequency*t
        ar=.025*self.frequency;uz=.03*self.frequency/k
        density=self.density*(1+.17*s+.1*(1-s)*np.cos(k*z)
            -(2*ar*(1-2*s)+k*uz*(1+.2*s))/self.frequency*np.cos(k*z)*np.sin(self.frequency*t))
        fr=ar*(1-s)*np.cos(k*z)*np.cos(self.frequency*t)
        fz=uz*(1+.2*s)*np.sin(k*z)*np.cos(self.frequency*t)+1e4*(1+.2*s)
        ft=.015*self.frequency*(1-s)*(1+.1*np.cos(phase))
        regular=np.stack((fr,ft,fz),axis=-1)*self.density/density[...,None]
        velocity=regular*np.stack((r,r,np.ones_like(r)),axis=-1)
        f=(1-s)**3;fs=-3*(1-s)**2;fss=6*(1-s)
        h=(1-s)**2;hs=-2*(1-s);primitive=s-s*s+s**3/3-1/3
        psi=.5*self.bias*r*r+self.poloidal*r*r*f*np.cos(phase)
        magnetic=np.stack((self.poloidal*k*r*f*np.sin(phase),r*self.toroidal*h*np.sin(other),
                           self.bias+2*self.poloidal*(f+s*fs)*np.cos(phase)),axis=-1)
        current=np.stack((-r*self.toroidal*k*h*np.cos(other),
            self.poloidal*r*(k*k*f-4/self.radius**2*(2*fs+s*fss))*np.cos(phase),
            2*self.toroidal*(h+s*hs)*np.sin(other)),axis=-1)/MU0
        electric=np.stack((np.zeros_like(r),-.7*self.frequency*self.poloidal*r*f*np.sin(phase),
            -.9*self.frequency*self.toroidal*self.radius**2/2*primitive*np.cos(other)),axis=-1)
        pe=self.density*E*self.electron_ev*(1+.13*s+.08*(1-s)*np.cos(k*z+.3*self.frequency*t))
        pi=self.density*E*self.ion_ev*(1+.11*s+.06*(1-s)*np.sin(k*z-.4*self.frequency*t))
        total=self.ion_mass+self.electron_mass;w=current/(density[...,None]*E)
        vi=velocity+self.electron_mass/total*w;ve=velocity-self.ion_mass/total*w
        return dict(n=density,v=velocity,regular_velocity=regular,j=current,b=magnetic,e=electric,
                    pe=pe,pi=pi,vi=vi,ve=ve,w=w,psi=psi,g=self.toroidal*h*np.sin(other))

    def derivatives(self,r,z,t):
        values=[]
        for axis,scale in enumerate((self.radius,self.length,1/self.frequency)):
            args=[r,z,t];step=1e-30*scale;args[axis]=np.asarray(args[axis],complex)+1j*step
            values.append({k:np.imag(v)/step for k,v in self.fields(*args).items()})
        return values

    def species_flux(self,r,z,t):
        a=self.fields(r,z,t);total=np.cross(a['e'],a['b'])/MU0
        for mass,velocity,pressure,gamma in ((self.ion_mass,a['vi'],a['pi'],self.gamma_i),
                                           (self.electron_mass,a['ve'],a['pe'],self.gamma_e)):
            total+=(.5*mass*a['n']*np.sum(velocity*velocity,axis=-1)+gamma*pressure/(gamma-1))[...,None]*velocity
        return total


def advection(a,b,dr,dz,r):
    result=a[...,0,None]*dr+a[...,2,None]*dz
    result[...,0]-=a[...,1]*b[...,1]/r
    result[...,1]+=a[...,1]*b[...,0]/r
    return result


def divergence(v,dr,dz,r):return dr[...,0]+v[...,0]/r+dz[...,2]


def parent_sources(ref,r,z,t):
    """Actual residual force/heat in each original species equation, no bulk formulas."""
    f=ref.fields(r,z,t);dr,dz,dt=ref.derivatives(r,z,t);forces={};heat={};continuity={}
    for key,mass,charge,gamma in (('i',ref.ion_mass,E,ref.gamma_i),('e',ref.electron_mass,-E,ref.gamma_e)):
        v=f['v'+key];p=f['p'+key];grad=np.stack((dr['p'+key],np.zeros_like(r),dz['p'+key]),axis=-1)
        acceleration=dt['v'+key]+advection(v,v,dr['v'+key],dz['v'+key],r)
        forces[key]=mass*f['n'][...,None]*acceleration-charge*f['n'][...,None]*(f['e']+np.cross(v,f['b']))+grad
        div=divergence(v,dr['v'+key],dz['v'+key],r)
        heat[key]=(dt['p'+key]+np.sum(v*grad,axis=-1)+gamma*p*div)/(gamma-1)
        continuity[key]=dt['n']+v[...,0]*dr['n']+v[...,2]*dz['n']+f['n']*div
    return forces,heat,continuity
