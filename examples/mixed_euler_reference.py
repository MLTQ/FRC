"""Time-dependent, boundary-compatible original cylindrical Euler manufacture."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E
from .log_fluid_reference import LogFluidReference


@dataclass(frozen=True)
class MixedEulerReference(LogFluidReference):
    kind: str='mixed_manufactured'

    @property
    def duration(self):return self.radius/self.sound

    def fields(self,r,z,time=0.):
        r,z=np.broadcast_arrays(r,z);y=(r/self.radius)**2;tau=self.sound*time/self.radius
        phase=np.pi*z/self.half_length-.4*tau;c=self.sound;n0=self.number_density
        n=n0*(1+.25*y)*(1+.2*np.cos(phase))
        omega=.2*c/self.radius*(1+.1*np.sin(phase))*(1+.05*np.cos(tau))
        vr=.15*c*r/self.radius*(1-y)*(1+.1*np.sin(phase))*(1+.1*np.cos(tau))
        vz=c*(.13+.2*(1+.15*y)*np.sin(phase));v=np.stack((vr,omega*r,vz),axis=-1)
        rotating=.5*self.mass*n0*omega**2*self.radius**2*(y+.125*y*y)*(1+.2*np.cos(phase))
        p=self.pressure*(1+.1*np.sin(phase))*(1+.2*(1-y)**2)*(1+.05*np.sin(tau))+rotating
        pe=p*(.6+.05*np.cos(phase+.2*np.sin(tau)));pi=p-pe
        se=np.log(pe/(n*E*(n/1e20)**(self.gamma_e-1)));si=np.log(pi/(n*E*(n/1e20)**(self.gamma_i-1)))
        return dict(n=n,velocity=v,electron=se,ion=si,pe=pe,pi=pi,Te=pe/(n*E),Ti=pi/(n*E),
                    energy=.5*self.mass*n*np.sum(v*v,axis=-1)+pe/(self.gamma_e-1)+pi/(self.gamma_i-1))

    def state(self,r,z,time=0.):
        f=self.fields(r,z,time);n=f['n'];u=np.concatenate((n[...,None],n[...,None]*f['velocity'],
                            n[...,None]*np.stack((f['electron'],f['ion']),axis=-1)),axis=-1)
        u[...,2]*=r
        return u

    def flux(self,r,z,time,direction):
        f=self.fields(r,z,time);flux=self.state(r,z,time)*f['velocity'][...,direction,None]
        flux[...,1+direction]+=(f['pe']+f['pi'])/self.mass
        return flux

    def source(self,r,z,time):
        """Continuous flux divergence, evaluated without production operators."""
        h=1e-28;f=self.fields(r,z,time)
        complex_r=r+1j*h;fr=self.fields(complex_r,z,time)
        dr=np.imag(complex_r[...,None]*self.state(complex_r,z,time)*fr['velocity'][...,0,None])/h
        dz=np.imag(self.flux(r,z+1j*h,time,2))/h
        dt=np.imag(self.state(r,z,time+1j*h))/h
        result=dt+dr/r[...,None]+dz
        # Differentiate pressure directly instead of subtracting p/r from
        # the derivative of r*p; the original equations are unchanged.
        result[...,1]+=np.imag(fr['pe']+fr['pi'])/(h*self.mass)-f['n']*f['velocity'][...,1]**2/r
        return result

    def primitive_source(self,r,z,time):
        """Separate primitive/material-equation derivation to check the source."""
        h=1e-28;f=self.fields(r,z,time)
        dr,dz,dt=({k:np.imag(v)/h for k,v in self.fields(rr,zz,tt).items()}
                  for rr,zz,tt in ((r+1j*h,z,time),(r,z+1j*h,time),(r,z,time+1j*h)))
        n=f['n'];v=f['velocity'];div=dr['velocity'][...,0]+v[...,0]/r+dz['velocity'][...,2]
        nt=-v[...,0]*dr['n']-v[...,2]*dz['n']-n*div
        vt=-v[...,0,None]*dr['velocity']-v[...,2,None]*dz['velocity']
        vt[...,0]+=v[...,1]**2/r-(dr['pe']+dr['pi'])/(self.mass*n)
        vt[...,1]-=v[...,0]*v[...,1]/r
        vt[...,2]-=(dz['pe']+dz['pi'])/(self.mass*n)
        result=np.concatenate(((dt['n']-nt)[...,None],(dt['n']-nt)[...,None]*v+n[...,None]*(dt['velocity']-vt),np.zeros(n.shape+(2,))),axis=-1)
        result[...,2]*=r
        for a,key,entropy,gamma in ((0,'pe','electron',self.gamma_e),(1,'pi','ion',self.gamma_i)):
            pt=-v[...,0]*dr[key]-v[...,2]*dz[key]-gamma*f[key]*div
            st=pt/f[key]-gamma*nt/n
            result[...,4+a]=(dt['n']-nt)*f[entropy]+n*(dt[entropy]-st)
        return result
