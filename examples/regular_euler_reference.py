"""Independent original Euler equations in nonsingular x=r²/2 coordinates."""
import numpy as np
from scipy.constants import elementary_charge as E
from .mixed_euler_reference import MixedEulerReference


class RegularEulerReference(MixedEulerReference):
    def fields_x(self,x,z,time=0.):
        x,z=np.broadcast_arrays(x,z);y=2*x/self.radius**2;tau=self.sound*time/self.radius;phase=np.pi*z/self.half_length-.4*tau;c=self.sound
        n=self.number_density*(1+.25*y)*(1+.2*np.cos(phase))
        u=.15*c/self.radius*(1-y)*(1+.1*np.sin(phase))*(1+.1*np.cos(tau))
        omega=.2*c/self.radius*(1+.1*np.sin(phase))*(1+.05*np.cos(tau));vz=c*(.13+.2*(1+.15*y)*np.sin(phase))
        rotating=.5*self.mass*self.number_density*omega**2*self.radius**2*(y+.125*y*y)*(1+.2*np.cos(phase))
        p=self.pressure*(1+.1*np.sin(phase))*(1+.2*(1-y)**2)*(1+.05*np.sin(tau))+rotating
        pe=p*(.6+.05*np.cos(phase+.2*np.sin(tau)));pi=p-pe
        entropy=np.log(np.stack((pe,pi),axis=-1)/(n[...,None]*E*(n[...,None]/1e20)**(np.array([self.gamma_e,self.gamma_i])-1)))
        return dict(n=n,u=u,omega=omega,vz=vz,pressure=np.stack((pe,pi),axis=-1),entropy=entropy)

    def regular_state(self,x,z,time=0.):
        f=self.fields_x(x,z,time);n=f['n']
        return np.concatenate((n[...,None],n[...,None]*np.stack((f['u'],f['omega'],f['vz']),axis=-1),n[...,None]*f['entropy']),axis=-1)

    def parent_rate(self,x,z,time=0.):
        f=self.fields_x(x,z,time);h=1e-28
        dx,dz=({k:np.imag(v)/h for k,v in self.fields_x(xx,zz,time).items()} for xx,zz in ((x+1j*h,z),(x,z+1j*h)))
        n=f['n'];u=f['u'];om=f['omega'];vz=f['vz'];speed=2*x*u
        div=2*u+2*x*dx['u']+dz['vz'];nt=-speed*dx['n']-vz*dz['n']-n*div
        ut=-speed*dx['u']-vz*dz['u']-u*u+om*om-np.sum(dx['pressure'],axis=-1)/(self.mass*n)
        ot=-speed*dx['omega']-vz*dz['omega']-2*u*om
        zt=-speed*dx['vz']-vz*dz['vz']-np.sum(dz['pressure'],axis=-1)/(self.mass*n)
        st=-speed[...,None]*dx['entropy']-vz[...,None]*dz['entropy']
        return np.concatenate((nt[...,None],np.stack((nt*u+n*ut,nt*om+n*ot,nt*vz+n*zt),axis=-1),nt[...,None]*f['entropy']+n[...,None]*st),axis=-1)

    def regular_source(self,x,z,time):
        return np.imag(self.regular_state(x,z,time+1j*1e-28))/1e-28-self.parent_rate(x,z,time)


def initial(driver,reference):
    state=reference.regular_state(driver.x,driver.z);state[0,0,:,1:3]=0.;return state


def energy(state,x,reference):
    """Independent physical energy, with the inactive axis slots irrelevant."""
    n=state[...,0];v=state[...,1:4]/n[...,None];s=state[...,4:]/n[...,None];g=np.array([reference.gamma_e,reference.gamma_i])
    p=n[...,None]*E*np.exp(s)*(n[...,None]/1e20)**(g-1)
    return reference.mass*n*(x*(v[...,0]**2+v[...,1]**2)+.5*v[...,2]**2)+np.sum(p/(g-1),axis=-1)
