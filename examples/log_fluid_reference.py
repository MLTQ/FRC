"""Independent translating contact, nonlinear simple wave and rotating equilibrium."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import atomic_mass,electron_mass,elementary_charge as E


@dataclass(frozen=True)
class LogFluidReference:
    kind: str='contact'
    radius: float=.035
    half_length: float=.12
    number_density: float=1e18
    ion_mass: float=39.948*atomic_mass
    pressure: float=.4
    gamma_e: float=5/3
    gamma_i: float=1.4
    direction: float=1.
    amplitude: float=.1

    @property
    def mass(self):return self.ion_mass+electron_mass
    @property
    def sound(self):return np.sqrt((.6*self.gamma_e+.4*self.gamma_i)*self.pressure/(self.mass*self.number_density))
    @property
    def duration(self):
        if self.kind=='contact':return 2*self.half_length/(.4*self.sound)
        if self.kind=='acoustic':return .5*self.half_length/self.sound
        if self.kind=='rotation':return .15*self.radius/self.sound
        raise ValueError('Known independent fluid reference required')

    def fields(self,r,z,time=0.):
        r,z=np.broadcast_arrays(r,z);k=np.pi/self.half_length;n0=self.number_density
        zero=np.zeros_like(r+z+time);vr=zero;vt=zero
        if self.kind=='contact':
            vz=zero+self.direction*.4*self.sound;phase=k*(z-vz*time)
            n=n0*(1+.3*np.cos(phase));pe=self.pressure*(.6+.15*np.sin(phase));pi=self.pressure-pe
        elif self.kind=='acoustic':
            if self.gamma_e!=self.gamma_i:raise ValueError('Simple wave uses equal stated species gammas')
            gamma=self.gamma_e;c0=self.sound;speed=self.amplitude*c0
            origin=z-self.direction*c0*time
            for _ in range(16):
                u=self.direction*speed*np.sin(k*origin)
                origin-=(origin+(self.direction*c0+.5*(gamma+1)*u)*time-z)/(1+.5*(gamma+1)*self.direction*speed*k*np.cos(k*origin)*time)
            vz=self.direction*speed*np.sin(k*origin)
            c=c0+.5*(gamma-1)*self.direction*vz
            n=n0*(c/c0)**(2/(gamma-1));pe=.6*self.pressure*(n/n0)**gamma;pi=.4*self.pressure*(n/n0)**gamma
        elif self.kind=='rotation':
            n=zero+n0;omega=.4*self.sound/self.radius;vt=omega*r;vz=zero+.15*self.sound
            total=self.pressure+.5*self.mass*n0*omega**2*r*r;pe=.6*total;pi=.4*total
        else:raise ValueError('Known independent fluid reference required')
        v=np.stack((vr,vt,vz),axis=-1)
        momentum=np.stack((n*vr,n*vz,n*r*vt),axis=-1)
        se=np.log(pe/(n*E*(n/1e20)**(self.gamma_e-1)))
        si=np.log(pi/(n*E*(n/1e20)**(self.gamma_i-1)))
        energy=.5*self.mass*n*np.sum(v*v,axis=-1)+pe/(self.gamma_e-1)+pi/(self.gamma_i-1)
        return dict(n=n,momentum=momentum,velocity=v,electron=se,ion=si,pe=pe,pi=pi,
                    Te=pe/(n*E),Ti=pi/(n*E),energy=energy)

    def derivatives(self,r,z,time=0.):
        h=1e-28
        return tuple({key:np.imag(value)/h for key,value in self.fields(rr,zz,tt).items()}
            for rr,zz,tt in ((r+1j*h,z,time),(r,z+1j*h,time),(r,z,time+1j*h)))


def average(mesh,function,order=16):
    q,w=np.polynomial.legendre.leggauss(order)
    r=mesh.r[:,None,None,None]+.5*mesh.dr*q[None,None,:,None]
    z=mesh.z[None,:,None,None]+.5*mesh.dz*q[None,None,None,:]
    r,z=np.broadcast_arrays(r,z);value=function(r,z)
    weight=.25*mesh.dr*r*w[None,None,:,None]*w[None,None,None,:]/mesh.radial_measure[:,None,None,None]
    return np.sum(value*weight[(...,)+(None,)*(value.ndim-4)],axis=(2,3))


def initial(driver,reference):
    mesh=driver.mesh;n=average(mesh,lambda r,z:reference.fields(r,z)['n'])
    momentum=average(mesh,lambda r,z:reference.fields(r,z)['momentum'])
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');f=reference.fields(r,z)
    return driver.initialize(n,momentum,f['electron'],f['ion'])
