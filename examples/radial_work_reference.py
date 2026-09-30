"""Independent initial Euler rates for a smooth fully two-dimensional rotating flow."""
import numpy as np
from scipy.constants import elementary_charge as E
from .log_fluid_reference import LogFluidReference


class RadialWorkReference(LogFluidReference):
    def fields(self,r,z,time=0.):
        r,z=np.broadcast_arrays(r,z);y=(r/self.radius)**2;phase=np.pi*z/self.half_length;c=self.sound
        n=self.number_density*(1+.25*y)*(1+.2*np.cos(phase));omega=.2*c/self.radius
        velocity=np.stack((.15*c*r/self.radius*(1-y)*(1+.1*np.sin(phase)),omega*r,
                           .2*c*(1+.15*y)*np.sin(phase)),axis=-1)
        rotating=.5*self.mass*self.number_density*omega**2*self.radius**2*(y+.125*y*y)*(1+.2*np.cos(phase))
        pressure=self.pressure*(1+.1*np.sin(phase))*(1+.2*(1-y)**2)+rotating
        pe=pressure*(.6+.05*np.cos(phase));pi=pressure-pe
        se=np.log(pe/(n*E*(n/1e20)**(self.gamma_e-1)));si=np.log(pi/(n*E*(n/1e20)**(self.gamma_i-1)))
        momentum=np.stack((n*velocity[...,0],n*velocity[...,2],n*r*velocity[...,1]),axis=-1)
        return dict(n=n,momentum=momentum,velocity=velocity,electron=se,ion=si,pe=pe,pi=pi,
                    Te=pe/(n*E),Ti=pi/(n*E),energy=.5*self.mass*n*np.sum(velocity**2,axis=-1)+pe/(self.gamma_e-1)+pi/(self.gamma_i-1))

    def rates(self,r,z):
        f=self.fields(r,z);h=1e-28
        dr={k:np.imag(v)/h for k,v in self.fields(r+1j*h,z).items()}
        dz={k:np.imag(v)/h for k,v in self.fields(r,z+1j*h).items()}
        v=f['velocity'];div=dr['velocity'][...,0]+v[...,0]/r+dz['velocity'][...,2]
        nt=-v[...,0]*dr['n']-v[...,2]*dz['n']-f['n']*div
        vt=-v[...,0,None]*dr['velocity']-v[...,2,None]*dz['velocity']
        vt[...,0]+=v[...,1]**2/r-(dr['pe']+dr['pi'])/(self.mass*f['n'])
        vt[...,1]-=v[...,0]*v[...,1]/r
        vt[...,2]-=(dz['pe']+dz['pi'])/(self.mass*f['n'])
        result=dict(n=nt,velocity=vt)
        for key,gamma in (('pe',self.gamma_e),('pi',self.gamma_i)):
            result[key]=-v[...,0]*dr[key]-v[...,2]*dz[key]-gamma*f[key]*div
        for key,short,gamma in (('electron','e',self.gamma_e),('ion','i',self.gamma_i)):
            result[key]=result['p'+short]/f['p'+short]-gamma*nt/f['n']
            result['T'+short]=f['T'+short]*(result['p'+short]/f['p'+short]-nt/f['n'])
        result['energy']=self.mass*(.5*nt*np.sum(v*v,axis=-1)+f['n']*np.sum(v*vt,axis=-1))+result['pe']/(self.gamma_e-1)+result['pi']/(self.gamma_i-1)
        return result

    def derivatives(self,r,z,time=0.):
        h=1e-28
        spatial=tuple({key:np.imag(value)/h for key,value in self.fields(rr,zz).items()}
                      for rr,zz in ((r+1j*h,z),(r,z+1j*h)))
        return (*spatial,self.rates(r,z))
